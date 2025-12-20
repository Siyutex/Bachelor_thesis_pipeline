# infer pseudotime for each cell in adata and necessary preprocessing
# input is an aggregated (batch corrected) h5ad file (can also have cnv annotated already)

import scanpy as sc
import helper_functions as hf
import os
import numpy as np
from scipy.sparse import issparse
import py_monocle as monocle
import warnings
import scipy
import pandas as pd
from sklearn.mixture import GaussianMixture


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


def compute_pseudotime_monocle(adata, root_cell_idx, layer):
    # leanr graph first (takes umap and a set of clusters)
    # then pseudotime (takes learn graph output, root cells, umap)
    
    # set internal adata
    internal_adata = adata.copy()

    # get umap
    sc.pp.pca(internal_adata, n_comps=50, svd_solver="arpack", layer=layer)
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


def return_smoothed_expression(adata, flavor, layer):
    # takes internal adata with only one matrix in adata.X
    # applies sliding window smoothing per gene along pseudotime
        # use numpy.convole with box kernel (or scipy.uniform_filter1d, which is more efficient)
        # mathematically, a convolution witha box kernel is equivalent to a moving average
        # for the numy function this would mean in an example for 1 gene:
            # function a = array of expression levels sorted by pseudotime index (f(x) = expression_level(pseudotime))
            # function b = array of size window_size with values 1/window_size
            # at the edges, where the kernel array stretches beyond the domain of the signal array, numpy pads the signal array with 0
                # the scipy function lets you decide wether to pad with 0, or mirror neighboring values or repeat the closest singal value (the latter 2 are more acurate for pseudotime in biological cells (expression does not just drop to 0 in the cell, that would be a processing artifact))   
    # returns new adata with changes applied

    # define window size as franction of cells (amount of discrete values in pseudotime space)
    if adata.n_obs < 1000:
        warnings.warn("Less than 1000 cells in dataset, smoothing may be locally biased")
    window_size = int(0.05 * adata.n_obs) # define window size as fraction of cells and make sure its an integer value
    vprint(f"Adata had {adata.n_obs} cells, using window size: {window_size}")

    # order cells by pseudotime
    order = np.argsort(adata.obs[f"{flavor}_pseudotime"])
    if layer in adata.layers.keys():
        X = adata.layers[layer][order, :]
    elif layer == "X":
        X = adata.X[order, :]

    # make sure X is a numpy array
    if not isinstance(X, np.ndarray):
        X = X.to_numpy()

    # apply smoothing
    smoothed = scipy.ndimage.uniform_filter1d(X, size=window_size, axis=0, mode="reflect")

    # restore original order (incase anything relies on it)
    inverse_order = np.argsort(order) # order is a list of indeces, first value that appears in the list is the index of the cell with lowest pseudotime; inverse order is a list of indeces, first value is the index at of the cell that was originally the first cell
    smoothed = smoothed[inverse_order] # apply transform to restore original indexing

    # apply smoothed matrix to new adata
    adata_new = adata.copy()
    if layer in adata_new.layers.keys():
        del adata_new.layers[layer]
        adata_new.layers[layer] = smoothed
    elif layer == "X":
        del adata_new.X
        adata_new.X = smoothed

    return adata_new

def filter_switches(adata, layer, flavor, seed: int = 42, bic_threshold: int = 10, mean_threshold: float = 0.5) -> sc.AnnData:
    """
    Filters genes for Boolean GRN inference by comparing 1-component vs 2-component GMMs.
    Updates adata.var with 'is_switch' and 'gmm_bics'.
    (genes where 1 component fits better are constants, and not important for network dynamics, so filtered out-
     genes with =>2 actual components all fit better to the 2 component model, they are switches)
    """
    # 1. Setup storage for results
    # We use a dictionary to store BIC values and a list for the switch boolean
    bic_results = {}
    switch_status = {gene: False for gene in adata.var_names}
    internal_adata = adata.copy()

    # Sort expression matrix by pseudotime (using a view/copy for calculation)
    pseudotime_col = f"{flavor}_pseudotime"
    if pseudotime_col not in internal_adata.obs:
        raise ValueError(f"{pseudotime_col} not found in adata.obs")
        
    adata_sorted = internal_adata[np.argsort(internal_adata.obs[pseudotime_col])].copy() # np.argsort returns the indices so it can be used inside []

    # Extract expression data
    matrix = adata_sorted.layers[layer] if layer else adata_sorted.X
    if scipy.sparse.issparse(matrix):
        matrix = matrix.toarray()
    
    data = pd.DataFrame(matrix, index=adata_sorted.obs_names, columns=adata_sorted.var_names)

    # 2. Iterate through genes
    i = 0
    for gene in data.columns:
        vprint(f"Curretnly working on gene: {i} of {len(data.columns)}")
        i+=1
        X = data[gene].values.reshape(-1, 1)
        
        # k=1 GMM
        gmm1 = GaussianMixture(n_components=1, random_state=seed, n_init=1).fit(X)
        bic1 = gmm1.bic(X)
        
        # k=2 GMM
        gmm2 = GaussianMixture(n_components=2, random_state=seed, n_init=1).fit(X)
        bic2 = gmm2.bic(X)
        
        # Store BICs in a dictionary format for this gene
        bic_results[gene] = {'k1': round(bic1, 2), 'k2': round(bic2, 2)}
        
        # Logic check for switch
        passed_bic = bic2 < (bic1 - bic_threshold)
        is_switch = False
        
        if passed_bic:
            means = np.sort(gmm2.means_.flatten())
            # Minimum delta check
            if abs(means[1] - means[0]) > min(means[1], means[0]) * mean_threshold:
                is_switch = True
        
        switch_status[gene] = is_switch

    # 3. Map results back to the ORIGINAL adata object
    # We use .map() to ensure the order matches adata.var_names exactly
    internal_adata.var['is_switch'] = internal_adata.var_names.map(switch_status)
    # save bic results in seperate columns (h5ad cannot handle dictionaries as values in columns)
    adata.var['bic_k1'] = adata.var_names.map({g: b['k1'] for g, b in bic_results.items()})
    adata.var['bic_k2'] = adata.var_names.map({g: b['k2'] for g, b in bic_results.items()})

    # 4. filter and return
    internal_adata = internal_adata[:, internal_adata.var['is_switch'] == True]
    return internal_adata





def main(input_data_file, output_data_dir, origin_clade, flavor, layer, smoothe_expression, find_switches, bic_threshold, mean_threshold):

    # define seed for reproducibility
    np.random.seed(42)

    # define vprint
    vprint = hf.make_vprint(verbose)


    adata = sc.read_h5ad(input_data_file)

    # import adata
    internal_adata = adata.copy()   

    # annotate root cell
    print("Annotating root cell...")
    root_idx = get_root_cell_clade(internal_adata, origin_clade)
    adata.uns["iroot"] = root_idx
    internal_adata.uns["iroot"] = root_idx

    # add pseudotime to internal adata
    print("Computing pseudotime...")
    if flavor == "monocle":
        compute_pseudotime_monocle(internal_adata, root_idx, layer)
        adata.obs["monocle_pseudotime"] = internal_adata.obs["monocle_pseudotime"]
    elif flavor == "dpt":
        compute_pseudotime_dpt(internal_adata)
        adata.obs["dpt_pseudotime"] = internal_adata.obs["dpt_pseudotime"]

    if smoothe_expression == True:
        print("Smoothing expression...")
        internal_adata = return_smoothed_expression(internal_adata, flavor, layer)
        if layer in adata.layers.keys():
            del adata.layers[layer]
            adata.layers[layer] = internal_adata.layers[layer]
        elif layer == "X":
            del adata.X
            adata.X = internal_adata.X

    if find_switches == True:
        print("Finding switch genes...")
        internal_adata = filter_switches(internal_adata, layer, flavor, seed=42, bic_threshold=bic_threshold, mean_threshold=mean_threshold)
        vprint(f"{internal_adata.shape[1]} switch genes found")
        adata = internal_adata.copy()

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))

if __name__ == "__main__":
    # import cmd args
    input_data_file, output_data_dir, origin_clade, flavor, layer, smoothe_expression, find_switches, bic_threshold, mean_threshold, verbose = hf.import_cmd_args(10)
    vprint = hf.make_vprint(verbose)

    main(input_data_file, output_data_dir, origin_clade, flavor, layer, smoothe_expression, find_switches, bic_threshold, mean_threshold)

    
