import scanpy as sc
import numpy as np

adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\CNV_inferred_PDAC.h5ad")

# print all keys and the contained datatypes
for key in adata.obs.keys():
    print(f"Obs annotation: {key} of datatype {type(adata.obs[key].iloc[0])}")
for key in adata.var.keys():
    print(f"Var annotation: {key} of datatype {type(adata.var[key].iloc[0])}")
for key in adata.obsm.keys():
    print(f"Obsm annotation: {key} of datatype {type(adata.obsm[key])} containing {type(adata.obsm[key][0,0])}")
for key in adata.uns.keys():
    print(f"Uns annotation: {key} of datatype {type(adata.uns[key])}")





