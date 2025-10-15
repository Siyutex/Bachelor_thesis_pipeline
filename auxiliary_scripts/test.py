import scanpy as sc

adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\cell_type_annotated\cell_type_annotated__PDAC_cancerous_0.h5ad")

# print all keys and the contained datatypes
for key in adata.obs.keys():
    print(f"Obs annotation: {key} of datatype {type(adata.obs[key][0])}")
for key in adata.var.keys():
    print(f"Var annotation: {key} of datatype {type(adata.var[key][0])}")


