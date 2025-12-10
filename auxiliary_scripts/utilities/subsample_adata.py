import scanpy as sc
import helper_functions as hf
import numpy as np
import os

def subsample_cells(file_path, layer, fraction, n_samples, output_dir):
    """
    Produce new adata with layer as X, then subsample to fraction of cells without replacement
    """

    # create isolated adata
    adata = sc.read(file_path)
    adata = hf.matrix_to_anndata(adata, layer).copy()

    # create boolean mask
    n_true = int(fraction * adata.shape[0]) # number of cells to keep
    mask = np.zeros(adata.shape[0], dtype=bool)
    mask[:n_true] = True # array with size n_obs, but ordered

    
    # subsample
    for i in range(n_samples):
        np.random.shuffle(mask) # shuffle to keep a random set of cells
        ss_adata = adata[mask,:].copy()
        ss_adata.write(os.path.join(output_dir, f"sample_{i}.h5ad"), compression="gzip")
        del ss_adata # free up memory
    

