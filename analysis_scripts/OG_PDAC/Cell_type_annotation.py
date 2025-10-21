# INPUT: preprocessed h5ad file 
# normalize the expression per cell (so score magnitude is not biased by overall expression magnitude of a cell)
# get z scores for each gene
# get putative cell annotation scores
# transfer cell type annotations to adata.obs

import scanpy as sc
import os 
import json
import numpy as np
from scipy import sparse
import pandas as pd
import helper_functions as hf
import scvi



def annotate_markers_z_score(adata, cutoff_unsure, cutoff_other): 
    # cutoff_unsure is a value between 0 and 1, specifying how high the second highest score can at most be relative to the highest to still annotate a well defined cell type
    # cutoff_other speciefies how many stdevs above / below the mean the lowest cell type score has to be to annotate a well defined cell type
    # normalize (if we don't, then cells with very high overall expression shift, the means of the genes
    # and thus the z score)
    # logarithmization not needed, because highly and lowly expressed markers already contribute on the same
    # scale because of 0 mean and counting stdevs instead of counts

    internal_adata = adata.copy() # so we don't modify the .X of the real data, because it is needed downstream
    
    if not hf.is_normalized(internal_adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(internal_adata, target_sum= 1e4)
    else:
        vprint("adata is already normalized")
    
    # get z scores for each gene
    vprint("Getting z scores for each gene...")
    X = np.array(internal_adata.X.todense())  # cells x genes, dense ndarray
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    internal_adata.X = sparse.csc_matrix((X-mu)/sigma)

    # get putative cell annotation scores (>0 means it likely is that cell type, <0 means it likely is not)
    vprint("Getting putative cell annotation scores...")
    for celltype, marker_list in MARKER_LISTS.items():
        if use_ensembl_ids:
            marker_levels = internal_adata[:, internal_adata.var["gene_symbols"].isin(marker_list)]
        else:
            marker_levels = internal_adata[:, internal_adata.var_names.isin(marker_list)]
        marker_levels = np.array(marker_levels.X.todense())

        # check if array is not empty (if it is empty that means that there is no overlap between
        #  marker gene set and the var names of the batch, ie that gene was removed in the batch
        #  because it was present in an insiginificant amount of cells)
        if marker_levels.size != 0:
            adata.obs[celltype + "_score"] = marker_levels.mean(axis=1) # want to add scores to actual adata
        else:
            adata.obs[celltype + "_score"] = np.zeros(adata.shape[0])

    # negative selection
    vprint("Correcting scores by negative markers...")
    for celltype, negative_marker_list in NEGATIVE_MARKER_LISTS.items():
        if use_ensembl_ids:
            negative_marker_levels = internal_adata[:, internal_adata.var["gene_symbols"].isin(negative_marker_list)]
        else:
            negative_marker_levels = internal_adata[:, internal_adata.var_names.isin(negative_marker_list)]
        negative_marker_levels = np.array(negative_marker_levels.X.todense())

        if marker_levels.size != 0:
            adata.obs[celltype + "_score"] -= marker_levels.mean(axis=1) # want to add scores to actual adata
        else:
            adata.obs[celltype + "_score"] -= np.zeros(adata.shape[0])

    adata.obs['cell_type'] = None
    score_columns = [key + "_score" for key in MARKER_LISTS.keys()]

    vprint("Transfering cell type annotations...")
    for cell in range(adata.shape[0]):
        ranked_score_list = []
        # rank cell type scores for the cell
        for column in score_columns:
            ranked_score_list.append((column.removesuffix("_score"), adata.obs.iloc[cell, adata.obs.columns.get_loc(column)])) # append a tuple with name of cell type and score for that cell type
        ranked_score_list = sorted(ranked_score_list, key=lambda x: x[1], reverse=True) # rank by the value at index 1 in the tuple (the actual score for that cell type)

        # Annotate
        # if all scores are below cutoff_other, annotate as "other"
        # if the second highest score is > cutoff_unsure * highest score, annotate as "unsure"
        # otherwise, annotate as the celltype with the highest score
        if any(score[1] > cutoff_other for score in ranked_score_list): 
            if ranked_score_list[1][1] > cutoff_unsure * ranked_score_list[0][1]:
                adata.obs.iloc[cell, adata.obs.columns.get_loc('cell_type')] = "unsure"
            else:
                adata.obs.iloc[cell, adata.obs.columns.get_loc('cell_type')] = ranked_score_list[0][0]
        else:
            adata.obs.iloc[cell, adata.obs.columns.get_loc('cell_type')] = "other"

    return None


def annotate_markers_cellassign(adata, use_ensembl_ids):
    
    vprint("Creating dataframe of marker lists...")
    # first make 1 pandas series per key (cell type) in the dict, list values (marker genes) will be the indexes, values will be 1
    marker_series_dict = {key: pd.Series(data=1, index=value) for key, value in MARKER_LISTS.items()}
    # now concatenate into dataframe (one column per series), fill missing values in each column (ie genes that are not markers for that cell type) with 0
    marker_df = pd.concat(marker_series_dict, axis=1).fillna(0)

    # create internal adata, to avoid modifying the real data
    internal_adata = adata.copy()

    # change internal adata var names to gene symbols, if ensembl ids were used (needed for cellassign)
    if use_ensembl_ids:
        internal_adata.var_names = internal_adata.var["gene_symbols"]

    # compute size factor (ratio of a cell's UMI counts to average cell's counts)
    lib_size = adata.X.sum(1)
    internal_adata.obs["size_factor"] = lib_size / np.mean(lib_size)

    # subset adata to only include marker genes
    vprint("Subsetting adata to only include marker genes...")
    internal_adata = internal_adata[:, internal_adata.var_names.isin(marker_df.index)].copy()
    vprint(f"Subset shape: {internal_adata.shape}")   

    # reformat X to csr matrix to make training faster
    internal_adata.X = sparse.csr_matrix(internal_adata.X)
   
    # assign a saved cellassign model if it exists
    try:
        model = scvi.external.CellAssign.load("my_cellassign_model/", internal_adata)
    except Exception as e:
        print("No fitting cellassign model found. Training a new one...")
        print(f"Model loading failed because of Exception: {e}")
        model = None
        pass

    # otherwise, train a new cellassign model
    if model is None:
        vprint("Setting up cellassign...")
        scvi.external.CellAssign.setup_anndata(internal_adata, size_factor_key="size_factor")
        model = scvi.external.CellAssign(internal_adata, cell_type_markers=marker_df)
        vprint("Training model...")
        model.train()
        model.save("my_cellassign_model/", overwrite=True)

    # predict cell types
    predictions = model.predict()

    # add predictions to adata
    adata.obs["cell_type"] = predictions.idxmax(axis=1).values
    print(adata.obs["cell_type"][:5])

    return None

def display_fractions(adata):
    # display global fraction of each cell type as dataframe
    data = []
    # go through your known marker-based cell types
    for cell_type, marker_list in MARKER_LISTS.items():
        frac = np.sum(adata.obs["cell_type"] == cell_type) / adata.shape[0]
        data.append((cell_type, frac))

    # add the "other" and "unsure" categories
    for special in ["other", "unsure"]:
        frac = np.sum(adata.obs["cell_type"] == special) / adata.shape[0]
        data.append((special, frac))

    # make DataFrame
    df = pd.DataFrame(data, columns=["cell_type", "fraction"]).round(3)
    print(df)

    if (df["fraction"] == 0).sum() >= df.shape[0]-1:
        raise ValueError("All but one cell types have a fraction of 0. Perhaps the use_ensembl_ids parameter does not match the input data?")


def save_outcome(adata):
    # save the processed data to temporary h5ad file, make relevant directory first
    final_output_path = os.path.join(output_data_dir, os.path.basename(input_data_file))
    adata.write(final_output_path, compression="gzip")

    # foward the path to the temporary file to the executor script via stdout
    print("Output: " + final_output_path)









if __name__ == "__main__":

    # import cmd args
    input_data_file, output_data_dir, use_ensembl_ids, marker_file_path, negative_marker_file_path, model, cutoff_unsure, cutoff_other, verbose = hf.import_cmd_args(9)
    vprint = hf.make_vprint(verbose)

    # set marker lists as global var
    with open(marker_file_path) as f:
        MARKER_LISTS = json.load(f)

    # set negative marker lists as global var, if present and needed
    if negative_marker_file_path is not None:
        with open(negative_marker_file_path) as f:
            NEGATIVE_MARKER_LISTS = json.load(f)
    elif negative_marker_file_path is None and model == "z_score":
        print("Running z_score without negative markers...")
        NEGATIVE_MARKER_LISTS = {}

    # import input data into anndata object (handles compressed h5ad files as well)
    print("Reading data...")
    adata = sc.read_h5ad(input_data_file)

    if model == "z_score":
        print("Annotating cell types using z score...")
        annotate_markers_z_score(adata, cutoff_unsure, cutoff_other) # second highest score has to be lower than cutoff unsure (eg 0.8 = 80%) of highest to be sure, any score has to be at least above cutoff other (eg 0.2) stdevs below the mean for that score among all cells
    elif model == "cellassign":
        print("Annotating cell types using Cellassign...")
        annotate_markers_cellassign(adata, use_ensembl_ids)

    print("Displaying fractions...")
    display_fractions(adata)
    
    print("Saving outcome...")
    save_outcome(adata)


