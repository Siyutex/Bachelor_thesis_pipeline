# This script isolates epithelial cells from a given mtx file

import scanpy as sc
import sys
import numpy as np

# import mtx file (usable as is in scanpy, industry standard)
input_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the mtx file.") #sys.argv takes command line arguments (or arguments from other scripts, here: biological_variance_pipeline_executor.py), first is the script name, second is the input file path
output_path = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to save the output file.")
adata = sc.read_h5ad(input_path)  # read preprocessed h5ad file, returns AnnData object


if 'EPCAM' in adata.var_names:  # check if EPCAM gene is present in the data
    epcam_expr = adata[:, 'EPCAM'].X # results in (cells, 1) column vector, where each cell has a value for EPCAM expression
    if hasattr(epcam_expr, "toarray"): 
        epcam_expr = epcam_expr.toarray().flatten() # flatten the array to a 1D array if it is a sparse matrix (which is the case for most scanpy data)
    else: 
        epcam_expr = np.array(epcam_expr).flatten() # flatten the array to a 1D array if it is not a sparse matrix

    threshold = 1.5  # cutoff chosen based on visual inspection of EPCAM expression in histogram & clustering of cells (done once to determine a good threshold)
    selected_cells = epcam_expr > threshold
    adata_epithelial = adata[selected_cells].copy()
    print(f"Number of epithelial cells isolated: {adata_epithelial.shape[0]}")
else:
    print(f"EPCAM gene not found in the data. No epithelial cells isolated from {input_path}.")
    adata_epithelial = adata.copy()