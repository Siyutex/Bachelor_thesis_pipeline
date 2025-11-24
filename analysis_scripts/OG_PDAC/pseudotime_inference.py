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
    # leanr graph first (takes umap and a set of clusters)
    # then pseudotime (takes learn graph output, root cells, umap)
    
    # set internal adata
    internal_adata = adata.copy()

    # get umap
    sc.pp.pca(internal_adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(internal_adata, n_neighbors=15, n_pcs=50, use_rep="X_pca")
    sc.tl.umap(internal_adata, min_dist=0.2) # same mindist as in plotting
    umap = internal_adata.obsm["X_umap"] # is a numpy ndarray

    # turn cnv clades into ndarray
    cnv_clades = internal_adata.obs["cnv_clade"].to_numpy()

    # learn principal graph
    projected_points, mst, centroids = monocle.learn_graph(matrix=umap, clusters=cnv_clades)
    
    # order cells along pseudotime
    pseudotime = monocle.order_cells(
        matrix=umap,
        centroids=centroids,
        mst=mst,
        projected_points=projected_points,
        root_cells=root_cell_idx
    )

    # add pseudotime to adata
    adata.obs["monocle_pseudotime"] = pseudotime



def main(input_data_file, output_data_dir, origin_clade, flavor):

    adata = sc.read_h5ad(input_data_file)

    # import adata
    internal_adata = hf.matrix_to_anndata(adata, "log1p")   

    # annotate root cell
    print("Annotating root cell...")
    root_idx = get_root_cell_clade(internal_adata, origin_clade)
    adata.uns["iroot"] = root_idx
    internal_adata.uns["iroot"] = root_idx

    # add pseudotime to internal adata
    if flavor == "monocle":
        compute_pseudotime_monocle(internal_adata, root_idx)
        adata.obs["monocle_pseudotime"] = internal_adata.obs["monocle_pseudotime"]
    elif flavor == "dpt":
        compute_pseudotime_dpt(internal_adata)
        adata.obs["dpt_pseudotime"] = internal_adata.obs["dpt_pseudotime"]

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))

if __name__ == "__main__":
    # import cmd args
    input_data_file, output_data_dir, origin_clade, flavor, verbose = hf.import_cmd_args(4)
    vprint = hf.make_vprint(verbose)

    main(input_data_file, output_data_dir, origin_clade, flavor)

    
