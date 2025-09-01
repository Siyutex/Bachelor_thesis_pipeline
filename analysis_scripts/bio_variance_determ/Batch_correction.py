# this scirpt is meant to be run as a subprocess of Biological_variance_pipeline_executor.py
# it will aggregate all datasets of a single tissue type (cancer and non cancer) into an anndata
# then annote all cells with lineage based on marker genes
# then run batch correction ber lineage (this will hopefully reduce batch effects aslo between cancer and non cancer cells)
# all of this will happen directly in gene expression space

import os
import sys
import tempfile
import scanpy as sc

# import command line arguments from ececutor script
input_data_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input directory")
output_data_path = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")


# generate anndata object from input h5ad files in the input directory
adata = None # initialize anndata object
for file in os.listdir(input_data_path):
    if file.endswith(".h5ad"):
        adata_new = sc.read_h5ad(os.path.join(input_data_path, file))
        
        # merge adata objects into one
        if adata is None:
            adata = adata_new
        else:
            adata = adata.concatenate(adata_new)