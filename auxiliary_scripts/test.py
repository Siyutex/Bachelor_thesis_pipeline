import scanpy as sc

adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\CNV_inferred_PDAC.h5ad")

X = adata.obsm("X_scANVI_corrected_gene_values_cnv").copy()

# check how many nans are in matrix
print(X.isna().sum().sum())

# check how many nans are in each row
print(X.isna().sum(axis=1))