# this script will produce 2D projections af the input file, either PCA or UMAP depending on the projection parameter

import os
import scanpy as sc
import matplotlib.pyplot as plt
import helper_functions as hf
import json
import sys
import infer_CNV as cnv
import numpy as np
import pandas as pd

"""def plot_with_external_legend(adata, projection: str, show=True, root_cell_idx: int = None, **kwargs):

    # Prevent Scanpy from immediately showing the plot  
    if projection == "UMAP":
        sc.tl.umap(adata)
        fig = sc.pl.umap(adata,  show=False, **kwargs)
    elif projection == "PCA":
        fig = sc.pl.pca(adata, show=False, **kwargs)
    ax = fig.axes

    if root_cell_idx != None:
        umap_coords = adata.obsm['X_umap']
        plt.scatter(umap_coords[root_cell_idx, 0], umap_coords[root_cell_idx, 1],
            color='red', s=100, edgecolor='black', label='highlighted cell')
        
        plt.show()
        return


    n_cols = 1
    # more than 15 legend entries -> 2 columns
    if len(ax.get_legend_handles_labels()[0]) > 20:
        n_cols = 2
    
    # Move legend outside and shrink font
    ax.legend(
        bbox_to_anchor=(1.05, 1),  # move right of plot
        loc='upper left',
        fontsize='small',
        title=kwargs.get('color', ''),  # use color key as legend title
        ncol=n_cols
    )
    
    # Adjust layout to avoid clipping
    plt.tight_layout()

    # always save to tempdir (which is output_data_dir)
    plt.savefig(os.path.join(output_data_dir, f"{projection}_colored_by_{kwargs.get('color', '')}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))

    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()"""


def matrix_to_anndata(adata, matrix_key):
    """createas a new anndata with the given matrix (any layer, obsm, or X) as adata.X and returns it.
       preserves obs annotations, only presever var annotatiosn if shape[1] matches
       
       DOES NOT MODIFY ADATA"""

    # check if matrix_key is unique among adata.layers, adata.obsm, and adata.X
    keys = list(adata.layers.keys()) + list(adata.obsm.keys()) + ["X"]
    if keys.count(matrix_key) > 1:
        raise ValueError(f"Matrix key '{matrix_key}' is not unique among layers, obsm, and X")

    # assign correct matrix
    if matrix_key    == "X":
        mat = adata.X.copy()
    elif matrix_key in adata.layers:
        mat = adata.layers[matrix_key].copy()
    elif matrix_key in adata.obsm:
        mat = adata.obsm[matrix_key].copy()
    else:
        raise ValueError(f"Key '{matrix_key}' not found in X, layers, or obsm")
    
    # assign obs and var
    obs = adata.obs.copy()
    var = adata.var.copy()

    # create new anndata, assign var if possible
    if len(var) == mat.shape[1]:
        vprint("Obsm matrix feature count matches adata.var, assigning var and obs to obsm_adata...")
        matrix_adata = sc.AnnData(X=mat, obs=obs, var=var)
    else:
        vprint("Obsm matrix feature count does not match adata.var, only assigning obs to obsm_adata...")
        matrix_adata = sc.AnnData(X=mat, obs=obs)

    return matrix_adata

def get_DEGs(adata, groupby):

    # rank genes groups overrides any previous contents of the fields it writes to 
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", use_raw=False, n_genes=10)
    
    # get text representation of top DEGs and save to tempdir
    rg = adata.uns["rank_genes_groups"]
    groups = rg["names"].dtype.names  # group names
    lines = []
    for group in groups:
        names = list(rg["names"][group])
        names = adata.var["gene_symbols"][names]
        names = ", ".join(names)
        lines.append(f"DEGs for {group}: {names}")
        if show:
            print(f"DEGs for {group}: {names}")
    with open(os.path.join(output_data_dir, f"DEGs_grouped_by_{groupby}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.txt"), "w") as f:
        f.write("\n".join(lines))

    # also plot as matrixplot and save to tempdir (CURRENTLY DOES NOT WORK, perhaps too many marker genes?)
    """if marker_file_path != None:
        vprint("Plotting DEG matrixplot...")
        sc.pl.matrixplot(adata, marker_genes_dict, groupby=groupby, dendrogram=False, cmap="viridis", standard_scale="var", colorbar_title="scaled \nexpression", show=False)
        plt.savefig(os.path.join(output_data_dir, f"DEGs_grouped_by_{groupby}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
        if show: 
            plt.show() # closes the plot, so run after savefig
            plt.close()"""


def plot_projection_and_DEGs(adata, layer):
    """
    layer should be something like "X" or "X_cnv" or "X_scvi_corrected"

    plots UMAP or PCA projection of most passed obs annotations (skips celltype if isolated) and leiden.
    also computes degs for cleaned obs annotations list.
    """
    # create temp anndata with layer as X, might not have var, if shape of layer and adata.Y don't match
    vprint(f"converting {layer} to anndata...")
    internal_adata = matrix_to_anndata(adata, layer)

    # normalize temp adata and store result in according field in internal adata 
    if not hf.is_normalized(internal_adata) and layer != "X_cnv": # X_cnv is log normalized, so this check does not work (special case)
        vprint("normalizing data...")
        sc.pp.normalize_total(internal_adata, target_sum=1e4)
    else:
        vprint(f"adata is already normalized, skipping to next step...")

    # pca on temp and store result in according field in internal adata
    vprint("computing PCA embeddings...")
    sc.pp.pca(internal_adata, svd_solver="arpack")

    vprint("finding neighbors...")
    sc.pp.neighbors(internal_adata, use_rep="X_pca", n_neighbors=15)

    vprint("computing leiden clusters...")
    sc.tl.leiden(internal_adata, resolution=0.1) # uses neighbor graph
    
    # internal list since lists are mutable and we do not want to change the global list
    colored_by = obs_annotations.copy() # smth like ["cell_type", "summed_cnv", "pseudotime", "cnv_score", ...]
    colored_by.append("leiden")
    
    # quality checks before plotting
    vprint("quality checks before plotting...")
    if "cell_type" in colored_by and cell_type is not None:
        colored_by.remove("cell_type")
        vprint(f"{cell_type} cells are isolated, removing \"cell_type\" from colored_by list")
    for entry in colored_by:
        if entry not in internal_adata.obs.keys():
            colored_by.remove(entry)
            vprint(f"{entry} is not in adata.obs, removing it from colored_by list")
    
    # determine suitable number of columns for figure (rougly 1.5 times as many columns as rows)
    n_cols = np.round(np.sqrt(len(colored_by)*1.5)).astype(int) 

    # create a figure with one plot per color, save to temp, show if show is True
    if projection == "UMAP":
        vprint("computing UMAP embedding...")
        sc.tl.umap(internal_adata) # uses neighbor graph
        sc.pl.umap(internal_adata, color=colored_by, show=False, ncols=n_cols)

        # highlight root cell
        if root_cell_idx != None:
            umap_coords = internal_adata.obsm['X_umap']
            plt.scatter(umap_coords[root_cell_idx, 0], umap_coords[root_cell_idx, 1],
                color='red', s=100, edgecolor='black', label='highlighted cell')

        plt.savefig(os.path.join(output_data_dir, f"UMAP_of_{layer}_({'all_celltypes' if cell_type == None else "just_" + cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
        
        if show: 
            plt.show()
        plt.close()
    elif projection == "PCA":
        sc.pl.pca(internal_adata, color=colored_by, show=False, ncols=n_cols)

        # highlight root cell
        if root_cell_idx != None:
            pca_coords = internal_adata.obsm['X_pca']
            plt.scatter(pca_coords[root_cell_idx, 0], pca_coords[root_cell_idx, 1],
                color='red', s=100, edgecolor='black', label='highlighted cell')

        plt.savefig(os.path.join(output_data_dir, f"PCA_of_{layer}_({'all_celltypes' if cell_type == None else "just_" + cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
        
        if show: 
            plt.show()
        plt.close()
    else:
        raise ValueError(f"projection {projection} is not supported")

    # for each grouping also run deg analysis
    if internal_adata.shape[1] == adata.shape[1]: # proxy for varnames being copied over (can't run DEG if you don't know which genes are present)
        
        vprint("Preparping data for DEG analysis...")
        # turn obs annotations into categorical, if they have a reasonable number of categories
        remove_colored_by = []
        for entry in colored_by:
            if len(internal_adata.obs[entry].unique()) <= 20:
                vprint(f"turning {entry} into categorical, unique values: {len(internal_adata.obs[entry].unique())}")
                internal_adata.obs[entry] = pd.Categorical(internal_adata.obs[entry])
            else:
                vprint(f"{entry} has too many categories ({len(internal_adata.obs[entry].unique())}), cannot compute DEGs for it")
                remove_colored_by.append(entry)

        colored_by = [x for x in colored_by if x not in remove_colored_by] # adjust list

        for entry in colored_by:
            vprint(f"computing DEGs for {entry}")
            get_DEGs(internal_adata, entry)
    else:
        vprint(f"skipping DEG analysis since {layer} has different number of genes than adata")


def pseudotime_vs_CNV(adata, show):

    sc.pl.scatter(adata, x="dpt_pseudotime", y="summed_cnvs", title="summed_cnv  vs pseudotime", show = False)
    # save the plot to tempdir
    plt.savefig(os.path.join(output_data_dir, f"summed_cnvs_vs_pseudotime_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
    # show the plot
    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()

    sc.pl.scatter(adata, x="dpt_pseudotime", y="cnv_score", title="cnv_score vs pseudotime", show = False)
    # save the plot to tempdir
    plt.savefig(os.path.join(output_data_dir, f"cnv_score_vs_pseudotime_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
    # show the plot
    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()



def main():

    # read data into anndata
    print("reading data...")
    adata = sc.read_h5ad(input_data_file)
    if verbose:
        vprint("adata summary:")
        vprint(adata) # show summary about adata structure
    
    # limit anndata to cell type, if not present, exit
    if cell_type is not None:
        print(f"isolating {cell_type} cells...")
        adata = adata[adata.obs["cell_type"] == cell_type].copy()
    if adata.X.shape[0] == 0:
        print(f"No {cell_type} found in data, skipping to next batch...")
        print("Output: " + "none")
        sys.exit(0)

    # UMAP / PCA plots + DEGs for each layer
    if "projections" in modules:
        for layer in layers: # layerlist like [adata.X, adata.obsm["X_cnv"], ...]
            if layer in adata.layers.keys() or layer in adata.obsm.keys() or layer == "X":
                print(f"plotting projections for {layer}...")
                plot_projection_and_DEGs(adata, layer=layer)
            else:
                raise ValueError(f"layer {layer} not found in adata.layers or adata.obsm or adata.X")
    else:
        vprint("skipping projections module...")

    # scatterplot of pseudotime vs summed_cnvs and cnv_score
    if "pseudotime_vs_cnv" in modules:
        print("plotting pseudotime_vs_cnv...")
        pseudotime_vs_CNV(adata, show=show)
    else:
        vprint("skipping pseudotime_vs_cnv module...")

    # OTHER MODULE IDEAS
    # force directed graph for pseudotime
    # DEG heatmap
    # phylogenetic tree for normal / transtion / tumor cells as per cnvscore
    # gene expression over pseudotime graph

    print(f"Output: {output_data_dir}")



if __name__ == "__main__":
    
    # import cmd args
    input_data_file, output_data_dir, modules, cell_type, obs_annotations, layers, projection, marker_file_path, root_cell_idx, show, verbose = hf.import_cmd_args(11)
    vprint = hf.make_vprint(verbose)

    if marker_file_path != None:
        with open(marker_file_path, "r") as f:
            marker_genes_dict = json.load(f)
    else:
        print("No marker file provided, skipping DEG matrix plots...")

    main()

    

