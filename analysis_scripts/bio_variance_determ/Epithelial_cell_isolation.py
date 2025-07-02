# This script isolates epithelial cells from a given mtx file

import scanpy
import sys

# import mtx file (usable as is in scanpy, industry standard)
input_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the mtx file.") #sys.argv takes command line arguments (or arguments from other scripts, here: biological_variance_pipeline_executor.py), first is the script name, second is the input file path
output_path = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to save the output file.")
adata = scanpy.read_mtx(input_path)  # read mtx file, returns AnnData object

# pca, neighbourhood graph construction, cluster
#scanpy.tl.pca(adata, n_comps=50)  # perform tl pca (tl = tookit, maintains original data, saves results in adata.obsm['X_pca']; pp = preprocessing, modifies adata.X)
#scanpy.pp.neighbors(adata, n_neighbors=10, n_pcs=50)  # construct neighbourhood graph
#scanpy.tl.leiden(adata, resolution=0.5)  # cluster with leiden algorithm

# map to 2D, color clusters according to marker genes
#scanpy.tl.umap(adata)  # compute UMAP embedding

# isolate epithelial cluster and put into intermediate file