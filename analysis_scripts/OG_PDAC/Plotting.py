# this script will produce 2D projections af the input file, either PCA or UMAP depending on the projection parameter

import os
import scanpy as sc
import matplotlib.pyplot as plt
import helper_functions as hf
import pandas as pd
import sys

def plot_with_external_legend(adata, projection: str, show=True, output_data_dir="", root_cell_idx: int = None, **kwargs):

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
    plt.close()

def get_DEGs(adata, groupby, cell_type, show):

    sc.tl.rank_genes_groups(adata, layer=obsm_layer, groupby=groupby, method="wilcoxon", use_raw=False, n_genes=10)
    
    rg = adata.uns["rank_genes_groups"]
    groups = rg["names"].dtype.names  # group names
    
    lines = []
    for group in groups:
        names = list(rg["names"][group])
        names = adata.var["gene_symbols"][names]
        names = ", ".join(names)
        lines.append(f"DEGs for {group}: {names}")
        vprint(f"DEGs for {group}: {names}")

    with open(os.path.join(output_data_dir, f"DEGs_grouped_by_{groupby}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.txt"), "w") as f:
        f.write("\n".join(lines))

    """# make a dataframe
    rg = adata.uns["rank_genes_groups"]
    groups = rg["names"].dtype.names  # group names

    data = {}
    for g in groups:
        names = list(rg["names"][g])
        if "ENS" in names[0]: # if the gene names are ensembl ids sawp them to gene symbols for readability
            names = adata.var["gene_symbols"][names]
        pvals = list(rg["pvals_adj"][g])
        data[g] = {"DEGs": names, "pvals_adj": pvals}

    df = pd.DataFrame.from_dict(data, orient="index")
    print(df)

    # plot the dataframe
    fig, ax = plt.subplots(figsize=(8, len(df) * 0.2))
    for i, (group, row) in enumerate(df.iterrows()):
        genes = ", ".join(row["DEGs"])
        ax.text(0, i, genes, fontsize=9, va="center")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df.index)
    ax.set_xlabel("Top DEGs")
    ax.set_title(f"Top DEGs grouped by {groupby}_(cell_type: {cell_type})")

    # save the plot to tempdir
    plt.savefig(os.path.join(output_data_dir, f"TOP_DEGs_grouped_by_{groupby}_(just_{cell_type})_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
    
    # show the plot
    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()"""


def pseudotime_vs_CNV(adata, show: bool = False, output_data_dir: str = ""):
    # spawn the plot
    sc.pl.scatter(adata, x="dpt_pseudotime", y="summed_cnvs", title="CNV vs pseudotime", show = False)

    # save the plot to tempdir
    plt.savefig(os.path.join(output_data_dir, f"CNV_vs_pseudotime_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))
    
    # show the plot
    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()



def main(input_data_file, output_data_dir, obs_annotations, obsm_layer, projection, cell_type, root_cell_idx, show):

    print("reading data...")
    adata = sc.read_h5ad(input_data_file)
    if obsm_layer == "X":
        adata.layers["X"] = adata.X
    else:
        adata.layers[obsm_layer] = adata.obsm[obsm_layer] # prbly doesn't work but worth a shot

    if cell_type is not None:
        print(f"isolating {cell_type} cells...")
        adata = adata[adata.obs["cell_type"] == cell_type].copy()

        if adata.X.shape[0] == 0:
            print(f"No {cell_type} found in data, skipping to next batch...")
            print("Output: " + "none")
            sys.exit(0)


    if not hf.is_normalized(adata, layer=obsm_layer):
        print("normalizing data...")
        sc.pp.normalize_total(adata, target_sum=1e4, layer=obsm_layer)
    else:
        print("Data already normalized, skipping to next step...")

    print("computing PCA embeddings...")
    sc.pp.pca(adata, svd_solver="arpack", layer=obsm_layer)

    print("finding neighbors...")
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)

    print("computing leiden clusters...")
    sc.tl.leiden(adata, resolution=0.1) # uses neighbor graph


    print(f"computing {projection} embeddings and plotting...")
    # one plot per annotation (eg cell type, cancer state, etc)
    for column in obs_annotations:
        get_DEGs(adata, column, cell_type, show=show)
        plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color=column, title=f"{projection} colored by {column} for {cell_type}")
    
    # also do leiden
    get_DEGs(adata, "leiden", cell_type, show=show)
    plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color="leiden", title=f"{projection} colored by leiden clusters for {cell_type}")

    # make pseudotime plot if pseudotime was computed
    if "dpt_pseudotime" in adata.obs.keys():
        pseudotime_vs_CNV(adata, show=show, output_data_dir=output_data_dir)
        plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color="dpt_pseudotime", color_map="plasma", vmin=0, vmax=1, title=f"{projection} colored by pseudotime for {cell_type}", root_cell_idx=root_cell_idx)
        plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color="summed_cnvs", color_map="viridis", title=f"{projection} colored by cnv counts for {cell_type}", root_cell_idx=root_cell_idx)

    print(f"Output: {output_data_dir}")



if __name__ == "__main__":
    
    # import cmd args
    input_data_file, output_data_dir, obs_annotations, obsm_layer, projection, cell_type, root_cell_idx, show, verbose = hf.import_cmd_args(7)
    vprint = hf.make_vprint(verbose)



    main(input_data_file, output_data_dir, obs_annotations, obsm_layer, projection, cell_type, root_cell_idx, show)

    

