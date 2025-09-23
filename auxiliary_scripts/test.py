import scanpy as sc
import os 
import sys
import numpy as np
import pandas as pd
from scipy.stats import zscore

adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\preprocessed\preprocessed_PDAC_cancerous_0.h5ad")

# global constants
THRESHOLD = 0.5
MARKER_LISTS = {
    "epithelial_cell": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "immune_cell": ["PTPRC"], # PTPRC = CD45
    "fibroblast": ["PDGFRA", "VIM", "COL1A1", "COL1A2", "LUM", "DCN"],
    "endothelial_cell": ["PECAM1", "VWF", "CLDN5"],
    "smooth_muscle_cell": ["ACTA2", "PDGFRB", "RGS5"],
}

THRESHOLD_Z = 1.0       # min z-score to call a cell type
REL_DOM_THRESHOLD = 0.3 # minimum relative dominance (top vs second best) to be confident
CTRL_SIZE = 50          # control gene size in score_genes


# normalize expression per cell
print("Normalizing expression")
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True) # does not affect n_counts in adata.obs, but affects adata.X (so even though each cell now has 10000 transcripts in X, n_counts is still what it was before normalization)
sc.pp.log1p(adata)




#-------------------------------------
# Annotation
#-------------------------------------


# ---------------------------
# 4. Compute scores
# ---------------------------
score_cols = []
for cell_type, genes in MARKER_LISTS.items():
    score_name = cell_type + "_score"
    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=score_name,
        ctrl_size=CTRL_SIZE
    )
    score_cols.append(score_name)

print(score_cols)

# ---------------------------
# 5. Create scores DataFrame
# ---------------------------
scores_df = adata.obs[score_cols]

print(scores_df.shape)
print(type(scores_df))

# ---------------------------
# 6. Per-cell z-score normalization across cell types
# ---------------------------
scores_z = scores_df.apply(zscore, axis=0)
print(type(scores_z))
print(scores_z.shape)
scores_z_df = pd.DataFrame(scores_z, index=scores_df.index, columns=scores_df.columns)

# ---------------------------
# 7. Determine top score and relative dominance
# ---------------------------
top_score = scores_z.max(axis=1)
second_score = scores_z.apply(lambda row: row.nlargest(2).iloc[-1], axis=1)
rel_dom = (top_score - second_score) / (top_score + 1e-6)  # relative difference
top_cell_type = scores_z.idxmax(axis=1).str.replace("_score", "")

# ---------------------------
# 8. Assign cell types
# ---------------------------
adata.obs["cell_type"] = np.where(
    top_score < THRESHOLD_Z,              # too low -> "other"
    "other",
    np.where(
        rel_dom < REL_DOM_THRESHOLD,      # top not dominant enough -> "unsure"
        "unsure",
        top_cell_type                      # confident assignment
    )
)

# debugging 
# get global lowest, highest and average score for each cell type
print("Global lowest, highest and average scores for each cell type")
for cell_type, marker_list in MARKER_LISTS.items():
    print(f"{cell_type}: lowest score: {np.min(adata.obs[cell_type + '_score'])}, highest score: {np.max(adata.obs[cell_type + '_score'])}, average score: {np.mean(adata.obs[cell_type + '_score'])}")

# get global relative amounts of each cell type 
print("Global relative amounts of each cell type")
for cell_type, marker_list in MARKER_LISTS.items():
    print(f"{cell_type}: {np.sum(adata.obs["cell_type"] == cell_type) / adata.shape[0]}")
print(f"Other: {np.sum(adata.obs['cell_type'] == 'other') / adata.shape[0]}")
