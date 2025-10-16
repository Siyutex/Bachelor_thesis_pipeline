# this script will produce 2D projections af the input file, either PCA or UMAP depending on the projection parameter

import os
import scanpy as sc
import matplotlib.pyplot as plt
import helper_functions as hf
import pandas as pd

def plot_with_external_legend(adata, projection: str, show=True, output_data_dir="", **kwargs):

    # Prevent Scanpy from immediately showing the plot  
    if projection == "UMAP":
        fig = sc.pl.umap(adata,  show=False, **kwargs)
    elif projection == "PCA":
        fig = sc.pl.pca(adata, show=False, **kwargs)
    ax = fig.axes

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
    plt.savefig(os.path.join(output_data_dir, f"{projection}_colored_by_{kwargs.get('color', '')}_for_{os.path.basename(input_data_file).removesuffix('.h5ad')}.png"))

    if show == True: plt.show() # closes the plot, so run after savefig
    plt.close()

def get_DEGs(adata, groupby):
    sc.tl.rank_genes_groups(adata, layer=layer, groupby=groupby, method="wilcoxon", use_raw=False, n_genes=10)
    df = pd.DataFrame(adata.uns["rank_genes_groups"]["pvals"])
    print(df)


def main(input_data_file, output_data_dir, obs_annotations, layer, projection, show):

    print("reading data...")
    adata = sc.read_h5ad(input_data_file)
    adata.layers["X"] = adata.X

    if not hf.is_normalized(adata, layer=layer):
        print("normalizing data...")
        sc.pp.normalize_total(adata, target_sum=1e4, layer=layer)
    else:
        print("Data already normalized, skipping to next step...")

    print("computing PCA embeddings...")
    sc.pp.pca(adata, svd_solver="arpack", layer=layer)

    print("finding neighbors...")
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)

    print("computing leiden clusters...")
    sc.tl.leiden(adata, resolution=0.1) # uses neighbor graph

    if projection == "UMAP":
        print("computing UMAP embeddings...")
        sc.tl.umap(adata) # uses neighbor graph

        # one plot per annotation (eg cell type, cancer state, etc)
        for column in obs_annotations:
            get_DEGs(adata, column)
            plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color=column, title=f"UMAP colored by {column} for {os.path.basename(input_data_file)}")
        
        # also do leiden
        get_DEGs(adata, "leiden")
        plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color="leiden", title=f"UMAP colored by leiden clusters for {os.path.basename(input_data_file)}")

    elif projection == "PCA":
        print("computing PCA embeddings...")

        # one plot per annotation (eg cell type, cancer state, etc)
        for column in obs_annotations:
            get_DEGs(adata, column)
            plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color=column, title=f"UMAP colored by {column} for {os.path.basename(input_data_file)}")
        
        # also do leiden
        get_DEGs(adata, "leiden")
        plot_with_external_legend(adata, projection, show=show, output_data_dir=output_data_dir, color="leiden", title=f"UMAP colored by leiden clusters for {os.path.basename(input_data_file)}")


    # plotting

    print(f"Output: {output_data_dir}")



if __name__ == "__main__":
    
    # import cmd args
    input_data_file, output_data_dir, obs_annotations, layer, projection, show, verbose = hf.import_cmd_args(7)
    vprint = hf.make_vprint(verbose)



    main(input_data_file, output_data_dir, obs_annotations, layer, projection, show)

    

