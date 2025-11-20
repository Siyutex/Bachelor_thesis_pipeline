# infer pseudotime for each cell in adata and necessary preprocessing
# input is an aggregated (batch corrected) h5ad file (can also have cnv annotated already)

import scanpy as sc
import helper_functions as hf
import os
import numpy as np
from scipy.sparse import issparse
import py_monocle as monocle


def check_normalize(adata):
    # normalize if not already
    if not hf.is_normalized(adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(adata, target_sum=1e4)
    else:
        vprint("adata is already normalized, skipping normalization...")


def annotate_root_cell_mincnv(adata):
    """
    Returns index of root cell. Root cell will be the cell with lowest cnv sum (based on gene_values_cnv)
    """

    # filter adata, use copy() to create new objects, not vies; assign to same name so old adata gets overridden and garbage collected
    vprint("Filtering adata for non_cancerous ductal cells")
    adata = adata[adata.obs["cell_type"] == "ductal_cell", :].copy()
    adata = adata[adata.obs["cancer_state"] == "non_cancerous", :].copy()

    if issparse(adata.X):
        X = np.array(adata.X)
    elif type(adata.X) == np.ndarray:
        X = adata.X

    # filter out NaN genes (inferCNV does not assign some genes a CNV value, if it skips over them with the sliding window)
    vprint("Filtering out NaN genes...")
    nan_columns = np.isnan(X).all(axis=0)
    X = X[:, ~nan_columns]

    # compute mean cnv for each gene
    vprint("Computing mean cnvs per gene...")
    median_cnvs_per_gene = np.median(X, axis=0)
    
    # compute sum of abs deviation of cnv from mean cnv for each cell
    vprint("Computing sum of abs deviation of cnv from mean cnv per cell...")
    sum_abs_dev = np.sum(np.abs(X - median_cnvs_per_gene), axis=1)

    best_cell_idx = np.argmin(sum_abs_dev)
    vprint(f"Best cell index: {best_cell_idx}")

    return best_cell_idx


def get_root_cell_clade(adata, origin_clade):
    """
    return root cell index. Root cell is a random cell from the given clade
    """
    # make sure clade column exists
    if "cnv_clade" not in adata.obs.columns:
        raise ValueError("adata.obs must have a column named 'cnv_clade'")

    # make sure target clade exists
    if origin_clade not in adata.obs["cnv_clade"].unique():
        raise ValueError(f"Clade {origin_clade} does not exist in adata.obs['cnv_clade']")

    # pick random cell where adata.obs["cnv_clade"] == origin_clade
    root_cell_idx = np.random.choice(np.where(adata.obs["cnv_clade"] == origin_clade)[0])
    vprint(f"Root cell index: {root_cell_idx}")

    return root_cell_idx


def compute_pseudotime_dpt(adata, n_pcs: int = 50, n_neighbors: int = 15):

    # compute PCA embedding (needed for neihbor graph), adds adata.obsm["X_pca"]
    vprint("Computing PCA embedding...")
    sc.pp.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    # compute neighbors, provides neighborhood graph which is used by diffmap
    vprint("Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca") # uses same n_neighbous as in UMAP plots so represenation is accurate

    # compute diffmap (automatically uses default fields created by neighbors)
    # random walks from each cell to all other cells to get distances
    vprint("Computing diffmap...")
    sc.tl.diffmap(adata, n_comps=n_neighbors)

    # compute pseudotime (uses diffusion distances to get pseudotime, and automatically uses default fields created by neighbors)
    # adds annotations to adata.obs["dpt_pseudotime"]
    print("Computing pseudotime...")
    sc.tl.dpt(adata, n_dcs=n_neighbors)


def compute_pseudotime_monocle(adata, root_cell_idx):
    monocle.compute_cell_states()
    monocle.differential_expression_genes()
    monocle.learn_graph()
    monocle.pseudotime()
    monocle.regression_analysis()

    return



def main(input_data_file, output_data_dir, origin_clade):

    adata = sc.read_h5ad(input_data_file)

    # import adata
    internal_adata = adata.copy()

    # check normalization, if not already done, normalize
    print("Checking normalization...")
    check_normalize(internal_adata)

    # annotate root cell
    print("Annotating root cell...")
    root_idx = get_root_cell_clade(internal_adata, origin_clade)
    adata.uns["iroot"] = root_idx
    internal_adata.uns["iroot"] = root_idx

    # prepare (pca, neighbors, diffmap) for pseudotime inference
    print("Adding necessary fields for pseudotime inference...")
    compute_pseudotime_dpt(internal_adata)

    # add pseudotime to actual adata
    adata.obs["pseudotime"] = internal_adata.obs["dpt_pseudotime"]

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))

if __name__ == "__main__":
    # import cmd args
    input_data_file, output_data_dir, origin_clade, verbose = hf.import_cmd_args(4)
    vprint = hf.make_vprint(verbose)

    main(input_data_file, output_data_dir, origin_clade)

    
