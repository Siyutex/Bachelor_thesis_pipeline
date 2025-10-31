import scanpy as sc
import helper_functions as hf

adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\CNV_inferred_PDAC_ductal_cell.h5ad")

new_adata = hf.matrix_to_anndata(adata, "X_scANVI_corrected")

print(f"Var names of old and new match: {adata.var_names == new_adata.var_names}")

# now show extract of both var_names
print(adata.var_names[:20])
print(new_adata.var_names[:20])

