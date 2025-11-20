import scanpy as sc

FILE_LOCATION = r"/home/julian/Bachelor_thesis_pipeline/Data/output_storage/isolated/isolated_transition_clades_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected.h5ad"

adata = sc.read_h5ad(FILE_LOCATION)
print(adata)
print(adata.X.shape)
print(adata.layers["log1p"].shape)