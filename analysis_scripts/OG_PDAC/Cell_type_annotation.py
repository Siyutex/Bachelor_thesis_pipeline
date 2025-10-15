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
from matplotlib import pyplot as plt
import helper_functions as hf




def annotate_markers_z_score(adata, cutoff_unsure, cutoff_other): 
    # cutoff_unsure is a value between 0 and 1, specifying how high the second highest score can at most be relative to the highest to still annotate a well defined cell type
    # cutoff_other speciefies how many stdevs above / below the mean the lowest cell type score has to be to annotate a well defined cell type
    # normalize (if we don't, then cells with very high overall expression shift, the means of the genes
    # and thus the z score)
    # logarithmization not needed, because highly and lowly expressed markers already contribute on the same
    # scale because of 0 mean and counting stdevs instead of counts

    internal_adata = adata.copy() # so we don't modify the .X of the real data, because it is needed downstream
    vprint("Normalizing adata...")
    sc.pp.normalize_total(internal_adata, target_sum= 1e4)
    
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
    input_data_file, output_data_dir, marker_file_path, cutoff_unsure, cutoff_other, use_ensembl_ids, verbose = hf.import_cmd_args(5)
    vprint = hf.make_vprint(verbose)

    with open(marker_file_path) as f:
        MARKER_LISTS = json.load(f)

    # import input data into anndata object (handles compressed h5ad files as well)
    print("Reading data...")
    adata = sc.read_h5ad(input_data_file)

    print("Annotating cell types...")
    annotate_markers_z_score(adata, cutoff_unsure, cutoff_other) # second highest score has to be lower than cutoff unsure (eg 0.8 = 80%) of highest to be sure, any score has to be at least above cutoff other (eg 0.2) stdevs below the mean for that score among all cells
    

    print("Displaying fractions...")
    display_fractions(adata)
    
    print("Saving outcome...")
    save_outcome(adata)


