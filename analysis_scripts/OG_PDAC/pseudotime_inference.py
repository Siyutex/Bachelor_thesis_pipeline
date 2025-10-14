# infer pseudotime for each cell in adata and necessary preprocessing
# input is an aggregated (batch corrected) h5ad file (can also have cnv annotated already)

import scanpy as sc
import helper_functions as hf
import os
import numpy as np

if __name__ == "__main__":
    # import cmd args
    input_data_file, output_data_dir, verbose = hf.import_cmd_args(3)
    vprint = hf.make_vprint(verbose)

    PREFIX = "cnv_"

    # import adata
    adata = sc.read_h5ad(input_data_file)

    # normalize if not already
    if not hf.is_normalized(adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(adata, target_sum=1e4)
    else:
        vprint("adata is already normalized, skipping normalization...")
    
    # compute PCA embedding (needed for neihbor graph), adds adata.obsm["X_pca"]
    vprint("Computing PCA embedding...")
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack")

    # compute neighbors, provides neighborhood graph which is used by diffmap
    vprint("Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep="X_pca")

    # compute diffmap (automatically uses default fields created by neighbors)
    # random walks from each cell to all other cells to get distances
    vprint("Computing diffmap...")
    sc.tl.diffmap(adata, n_comps=15)

    # annotate root cells (do least cnvs, ductal, non cancerous and maybe some other criteria)
    # should be cancaer_state == non_cancerous, ductal_cell, and cnvs of choosen cell for all genes as close as possible to median cnvs of ductal+nc
    # CURRENT ISSUE: to many nans in gene_values_cnv bcs inferCNV needs densly packed gene set per chromosome (but we only use HVGs atm)
    vprint("Annotating root cell...")
    adata.uns["iroot"] = np.flatnonzero(adata.obs["cell_type"] == "ductal_cell")[0] # choose first ductal cell


    # fix smth like this
    """
    # filter adata, use copy() to create new objects, not vies; assign to same name so old adata gets overridden and garbage collected
adata = adata[adata.obs["cell_type"] == "ductal_cell", :].copy()
adata = adata[adata.obs["cancer_state"] == "non_cancerous", :].copy()

if issparse(adata.layers["gene_values_cnv"]):
    X = np.array(adata.layers["gene_values_cnv"])
elif type(adata.layers["gene_values_cnv"]) == np.ndarray:
    X = adata.layers["gene_values_cnv"]

# compute median cnv number for each gene in those cells
median_cnvs = np.median(X, axis=0)
print(f"Median cnvs: {median_cnvs[:5]}")

# now find cell where sum(abs(median_cnvs - cnvs)) is smallest for all genes
# compute per-cell L1 deviation from median

print((X - median_cnvs)[:5, :5])
print(np.abs(X - median_cnvs)[:5, :5])

deviations = np.sum(np.abs(X - median_cnvs), axis=1)
print(f"Deviations: {deviations[:5]}")

# find index of cell with smallest total deviation
best_cell_idx = np.argmin(deviations)
    """

    # compute pseudotime (uses diffusion distances to get pseudotime, and automatically uses default fields created by neighbors)
    # adds annotations to adata.obs["dpt_pseudotime"]
    vprint("Computing pseudotime...")
    sc.tl.dpt(adata, n_dcs=15)


    # TESTING
    sc.pl.scatter(adata, x="dpt_pseudotime", y="summed_cnvs", title="CNV vs pseudotime")


    # export adata
    vprint("Exporting adata...")
    adata.write(os.path.join(output_data_dir, "pseudotime_" + os.path.basename(input_data_file).removeprefix(PREFIX)), compression="gzip")

    # foward the path to the temporary file to the executor script via stdout
    print("Output: " + os.path.join(output_data_dir, "pseudotime_" + os.path.basename(input_data_file).removeprefix(PREFIX)))