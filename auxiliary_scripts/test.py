import scanpy as sc

adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\reduced\reduced_PDAC_ductal_cell.h5ad")

adata = adata[:200,:] # first 200 cells

adata.write(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\tree\test trees\partial_real_data\test.h5ad")