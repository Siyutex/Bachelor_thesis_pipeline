import scanpy as sc
import numpy as np
import helper_functions as hf
import os 
import pandas as pd

adata1 = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\CNV_inferred_PDAC_ductal_cell.h5ad")
adata2 = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\reduced\reduced_PDAC_ductal_cell_X_is_X_scANVI_corrected.h5ad")

# check numeric identiy
if np.allclose(adata1.obsm["X_scANVI_corrected"], adata2.X):
    print("Numeric identity check passed")
else:
    print("Numeric identity check failed")

# check type identity
if type(adata1.obsm["X_scANVI_corrected"]) == type(adata2.X):
    print("Type identity check passed")
else:
    print("Type identity check failed")
    print(f"Type of adata1: {type(adata1.obsm['X_scANVI_corrected'])}")
    print(f"Type of adata2: {type(adata2.X)}")

# check mean identtity 
if np.allclose(adata1.obsm["X_scANVI_corrected"].mean().mean(), adata2.X.mean()):
    print("Mean identity check passed")
else:
    print("Mean identity check failed")
    print(f"Mean of adata1: {adata1.obsm['X_scANVI_corrected'].mean()}")
    print(f"Mean of adata2: {adata2.X.mean()}")

# check order identity (cells in same order)
if np.all(adata1.obs_names == adata2.obs_names):
    print("Order identity check passed")
else:
    print("Order identity check failed")