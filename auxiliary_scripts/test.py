import scanpy as sc
import numpy as np
from scipy import sparse
import pandas as pd
import os


# just ductal vs acinar
MARKER_LISTS = { # verified with cellmarker 2.0
    # epithelial 
    "ductal_cell": ["CFTR", "SOX9", "MUC1", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "acinar_cell": ["PRSS1", "PRSS2", "CPA1", "CPA2", "AMY2A/B", "CELA3A", "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],

    # endocrine
    "α-cell": ["GCG"],
    "β-cell": ["INS"],
    "δ-cell": ["SST"],
    "PP-cell": ["PPY"],
    "ε-cell": ["GHRL"],

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
    # normlaize (if we don't, then cells with very high overall expression shift, the means of the genes
    # and thus the z score)
    # logarithmization not needed, because highly and lowly expressed markers already contribute on the same
    # scale because of 0 mean and counting stdevs instead of counts

    internal_adata = adata.copy() # so we don't modify the .X of the real data, because it is needed downstream
    sc.pp.normalize_total(internal_adata, target_sum= 1e4)
    
    # get z scores for each gene
    print(f"adata shape {internal_adata.shape}")
    X = np.array(internal_adata.X.todense())  # cells x genes, dense ndarray
    print(type(X))
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    internal_adata.X = sparse.csc_matrix((X-mu)/sigma)

    # get putative cell annotation scores (>0 means it likely is that cell type, <0 means it likely is not)
    for celltype, marker_list in MARKER_LISTS.items():
        marker_levels = internal_adata[:, internal_adata.var_names.isin(marker_list)]
        print("marker levels adata shape", marker_levels.shape)
        marker_levels = np.array(marker_levels.X.todense())
        print("marker levels array size", marker_levels.size)

        # check if array is not empty (if it is empty that means that there is no overlap between
        #  marker gene set and the var names of the batch, ie that gene was removed in the batch
        #  because it was present in an insiginificant amount of cells)
        if marker_levels.size != 0:
            adata.obs[celltype + "_score"] = marker_levels.mean(axis=1) # want to add scores to actual adata
        else:
            print("ELSE TIRGGERED")
            adata.obs[celltype + "_score"] = np.zeros(adata.shape[0])

    adata.obs['cell_type'] = None
    score_columns = [key + "_score" for key in MARKER_LISTS.keys()]

    for cell in range(adata.shape[0]):
        ranked_score_list = []
        for column in score_columns:
            ranked_score_list.append((column.removesuffix("_score"), adata.obs.iloc[cell, adata.obs.columns.get_loc(column)])) # append a tuple with name of cell type and score for that cell type
        ranked_score_list = sorted(ranked_score_list, key=lambda x: x[1], reverse=True) # rank by the value at index 1 in the tuple (the actual score for that cell type)

        if any(score[1] > cutoff_other for score in ranked_score_list): 
            if ranked_score_list[1][1] > cutoff_unsure * ranked_score_list[0][1]:
                adata.obs.iloc[cell, adata.obs.columns.get_loc('cell_type')] = "unsure"
            else:
                adata.obs.iloc[cell, adata.obs.columns.get_loc('cell_type')] = ranked_score_list[0][0]
        else:
            adata.obs.iloc[cell, adata.obs.columns.get_loc('cell_type')] = "other"

    return None




data_dir = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\preprocessed"


ductal_count = 0
acinar_count = 0
cell_count = 0
for file in os.listdir(data_dir):
    print("Reading data")
    adata = sc.read_h5ad(os.path.join(data_dir, file))

    annotate_markers_z_score(adata, 0.8, -0.2)

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


    # amount of epithelial cells
    print(f"Acinar cells: {np.sum(adata.obs['cell_type'] == 'acinar_cell')}")
    print(f"Ductal cells: {np.sum(adata.obs['cell_type'] == 'ductal_cell')}")

    ductal_count += np.sum(adata.obs['cell_type'] == 'ductal_cell')
    acinar_count += np.sum(adata.obs['cell_type'] == 'acinar_cell')
    cell_count += adata.shape[0]

print(f"Total ductal cells: {ductal_count}")
print(f"total ductal fraction: {ductal_count / cell_count}")
print(f"Total acinar cells: {acinar_count}")
print(f"total acinar fraction: {acinar_count / cell_count}")
print(f"Total immune cells: {np.sum(adata.obs['cell_type'] == 'immune_cell')}")
print(f"Total cells: {cell_count}")