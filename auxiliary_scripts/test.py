# get root cell for pseudotime
# should be cancaer_state == non_cancerous,
# ductal_cell, and cnvs as close as possible to avg cvns of ductal cells

import scanpy as sc
import numpy as np
from scipy.sparse import issparse
import json


adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\cnv_annotated_PDAC_cancerous_1.h5ad")
adata = adata[:100,:].copy()
adata.layers["gene_values_cnv"] = np.array(adata.layers["gene_values_cnv"])


print(f"type of gene_value counts is {type(adata.layers["gene_values_cnv"])}, shape {adata.layers["gene_values_cnv"].shape}")
nan_count = np.isnan(adata.layers["gene_values_cnv"]).sum()
print(f"Number of NaNs in gene_value counts: {nan_count}")
nan_columns = np.any(np.isnan(adata.layers["gene_values_cnv"]), axis=0) # true if any element in column is NaN
print(f"Columns with NaNs: {nan_columns}") 

nan_col_indices = np.where(nan_columns)[0]
with open("nan_col_indices.json", "w") as f:
    json.dump(nan_col_indices.tolist(), f)


print(f"type of X is {type(adata.X)}, shape {adata.X.shape}")
nan_count = np.isnan(adata.X).sum()
print(f"Number of NaNs in adata.X: {nan_count}")  

# filter adata, use copy() to create new objects, not vies; assign to same name so old adata gets overridden and garbage collected
adata = adata[adata.obs["cell_type"] == "ductal_cell", :].copy()
adata = adata[adata.obs["cancer_state"] == "non_cancerous", :].copy()

if issparse(adata.layers["gene_values_cnv"]):
    X = np.array(adata.layers["gene_values_cnv"])
elif type(adata.layers["gene_values_cnv"]) == np.ndarray:
    X = adata.layers["gene_values_cnv"]

# compute median cnv number for each gene in those cells
median_cnvs = np.median(X, axis=0)
print(f"Median cnvs: {median_cnvs[:5]}")

# now find cell where sum(abs(median_cnvs - cnvs)) is smallest for all genes
# compute per-cell L1 deviation from median

print((X - median_cnvs)[:5, :5])
print(np.abs(X - median_cnvs)[:5, :5])

deviations = np.sum(np.abs(X - median_cnvs), axis=1)
print(f"Deviations: {deviations[:5]}")

# find index of cell with smallest total deviation
best_cell_idx = np.argmin(deviations)

print(best_cell_idx)
print(adata.obs_names[best_cell_idx])
