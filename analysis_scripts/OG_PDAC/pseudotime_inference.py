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
    vprint("Annotating root cell...")
    adata.uns["iroot"] = np.flatnonzero(adata.obs["cell_type"] == "ductal_cell")[0] # choose first ductal cell

    # compute pseudotime (uses diffusion distances to get pseudotime, and automatically uses default fields created by neighbors)
    # adds annotations to adata.obs["dpt_pseudotime"]
    vprint("Computing pseudotime...")
    print(adata.uns.keys())
    sc.tl.dpt(adata, n_dcs=15)


    # TESTING
    sc.pl.scatter(adata, x="dpt_pseudotime", y="summed_cnvs", title="CNV vs pseudotime")


    # export adata
    vprint("Exporting adata...")
    adata.write(os.path.join(output_data_dir, "pseudotime_" + os.path.basename(input_data_file).removeprefix(PREFIX)), compression="gzip")

    # foward the path to the temporary file to the executor script via stdout
    print("Output: " + os.path.join(output_data_dir, "pseudotime_" + os.path.basename(input_data_file).removeprefix(PREFIX)))