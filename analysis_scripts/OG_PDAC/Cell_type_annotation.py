# INPUT: preprocessed h5ad file 
# normalize the expression per cell (so score magnitude is not biased by overall expression magnitude of a cell)
# log1p, so very highly expressed markers don't drown out signal of lowly / normally expressed markers
# sc.tl.score_genes for each set of markers, then choose cell type with highest score for each cell
# if one score is significantly higher than the others, choose that cell type
# if that is not the case choose "other"

import scanpy as sc
import os 
import sys
import numpy as np
from scipy import sparse
import pandas as pd
from matplotlib import pyplot as plt

# import command line arguments from ececutor script
input_data_file = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input file")
output_data_dir = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")

# import input data into anndata object (handles compressed h5ad files as well)
print("Reading data")
adata = sc.read_h5ad(input_data_file)

# original marker list
MARKER_LISTS = { # verified with cellmarker 2.0
    # epithelial 
    "ductal_cell": ["CFTR", "SOX9", "MUC1", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "acinar_cell": ["PRSS1", "PRSS2", "CPA1", "CPA2", "AMY2A/B", "CELA3A", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],

    # endocrine
    "alpha-cell": ["GCG"],
    "beta-cell": ["INS"],
    "delta-cell": ["SST"],
    "PP-cell": ["PPY"],
    "epsilon-cell": ["GHRL"],

    # immunte
    "T_cell": ["PTPRC", "CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell": ["PTPRC", "CD19", "CD79A", "MS4A1", "MZB1"],
    "NK_cell": ["PTPRC", "NCAM1", "NKG7", "GNLY"],
    "macrophage": ["PTPRC", "CD14", "CD68", "LYZ", "FCGR3A"],
    "dendritic_cell": ["PTPRC", "ITGAX", "CLEC9A"],
    "neutrophil": ["PTPRC", "S100A8", "S100A9", "FCGR3B"],
    "mast_cell": ["PTPRC", "TPSAB1", "KIT"],

    # stromal
    "fibroblast": ["PDGFRA", "VIM", "COL1A1", "COL1A2", "LUM", "DCN"],
    "endothelial_cell": ["PECAM1", "VWF", "CLDN5"],
    "smooth_muscle_cell": ["ACTA2", "PDGFRB", "RGS5"],
}

def annotate_markers_z_score(adata, cutoff_unsure, cutoff_other): # cutoff is a value between 0 and 1, specifying how high the second highest score can at most be relative to the highest to still annotate a well defined cell type
    # normalize (if we don't, then cells with very high overall expression shift, the means of the genes
    # and thus the z score)
    # logarithmization not needed, because highly and lowly expressed markers already contribute on the same
    # scale because of 0 mean and counting stdevs instead of counts

    internal_adata = adata.copy() # so we don't modify the .X of the real data, because it is needed downstream
    sc.pp.normalize_total(internal_adata, target_sum= 1e4)
    
    # get z scores for each gene
    X = np.array(internal_adata.X.todense())  # cells x genes, dense ndarray
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    internal_adata.X = sparse.csc_matrix((X-mu)/sigma)

    # get putative cell annotation scores (>0 means it likely is that cell type, <0 means it likely is not)
    for celltype, marker_list in MARKER_LISTS.items():
        marker_levels = internal_adata[:, internal_adata.var["gene_symbols"].isin(marker_list)]
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


def save_outcome(adata):
    # save the processed data to temporary h5ad file, make relevant directory first
    final_output_path = os.path.join(output_data_dir,  "annotated_" + os.path.basename(input_data_file).removeprefix("preprocessed_"))
    adata.write(final_output_path, compression="gzip")

    # foward the path to the temporary file to the executor script via stdout
    print("Output: " + final_output_path)









annotate_markers_z_score(adata, 0.8, -0.2) # second highest score has to be lower than 80% of highest to be sure, any score has to be at least above 0.2 stdevs below the mean for that score among all cells
display_fractions(adata)
save_outcome(adata)


