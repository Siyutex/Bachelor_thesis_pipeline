import scanpy as sc

FILE_LOCATION = r"/home/julian/Bachelor_thesis_pipeline/Data/output_storage/isolated/isolated_transition_clades_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected+log1p.h5ad"

adata = sc.read_h5ad(FILE_LOCATION)
print(f"Adata summary: \n{adata}")
print(f"Head of adata: \n{adata.X[:5]}")
print(f"library size: \n{adata.X.sum(axis=1)}")
for layer in adata.layers.keys():
    print(f"Layer {layer} summary: \n{adata.layers[layer]}")
