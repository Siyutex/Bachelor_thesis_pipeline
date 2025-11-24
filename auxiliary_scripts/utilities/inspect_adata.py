import scanpy as sc

FILE_LOCATION = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/pseudotime/pseudotime_inferred_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected_cancer_state_inferred_tree_is_['transitional'].h5ad"
layer_to_check = "X_scANVI_corrected"


adata = sc.read_h5ad(FILE_LOCATION)
print(f"Adata summary: \n{adata}")
if layer_to_check in adata.layers.keys():
    print(f"Library size of {layer_to_check} is {adata.layers[layer_to_check].sum(axis=1)}")
elif layer_to_check in adata.obsm.keys():
    print(f"Library size of {layer_to_check} is {adata.obsm[layer_to_check].sum(axis=1)}")
if adata.X != None:
    print(f"Head of adata.X: \n{adata.X[:5]}")
for layer in adata.layers.keys():
    print(f"Layer {layer} summary: \n{adata.layers[layer]}")
for layer in adata.obsm.keys():
    print(f"Obsm {layer} summary: \n{adata.obsm[layer]}")

