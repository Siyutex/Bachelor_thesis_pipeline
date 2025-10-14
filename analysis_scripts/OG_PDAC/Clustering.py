# this script will produce 2D embeddings with UMAP for visualization
# it can take any number of annotated batches as input
# if several batches are passed, it will merge them into 1 object, then cluster and plot, batch, leiden and annotations
# if one batch is passed, it will just cluster and plot, leiden and annotations 

import os
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import helper_functions as hf

# import cmd args
input_data_file_or_dir, output_data_dir, annotation_columns, embedding, verbose = hf.import_cmd_args(5)
vprint = hf.make_vprint(verbose)

def aggregate_batches(input_dir):
    """
    Create an anndata object from a list of h5ad files. Annotate batch origin, cancer state
    and keep any previous annotations.
    """

    batches = []
    batch_labels = []
    cancer_states = []

    vprint("Reading data")
    for file in os.listdir(input_dir):
            adata = sc.read_h5ad(os.path.join(input_dir, file))
            batches.append(adata)
            batch_labels.append(file) 
            cancer_states.append("non_cancerous" if "non_cancerous" in file else "cancerous")
        
    vprint("Cancer states of batches:", cancer_states)

    # Expand cancer states to match the number of cells in each batch
    vprint("Expanding cancer states to fit with adata.obs.shape[0]")
    cancer_states_expanded = []
    for state, ad in zip(cancer_states, batches):
        cancer_states_expanded.extend([state] * ad.n_obs)

    # Merge all batches into one AnnData
    vprint("Merging batches")
    merged = anndata.concat(
        batches,
        join="outer",           # union of genes, fill missing with 0
        label="batch",        # name of the new obs column
        keys=batch_labels,     # values for the new obs column
        fill_value=0,
        index_unique="-"
    )

    # Add cancer state information to .obs
    merged.obs["cancer_state"] = cancer_states_expanded

    # show structure of resulting adata object
    vprint(f"Structure of merged adata object: {merged}")

    return merged


def plot_with_external_legend(adata, show=True, save=False, **kwargs):
    # Prevent Scanpy from immediately showing the plot
    fig, ax = plt.subplots(figsize=(6,6))  
    sc.pl.umap(adata, ax=ax,  show=False, **kwargs)

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
    if show == True: plt.show()
    if save == True: plt.savefig(os.path.join(output_data_dir, f"UMAP colored by {kwargs.get('color', '')} for {os.path.basename(input_data_file_or_dir).removesuffix('.h5ad')}.png"))


# check if input is directory of file
if os.path.isdir(input_data_file_or_dir):
    vprint(f"Input is a directory with {len(os.listdir(input_data_file_or_dir))} files")
    vprint("aggregating batches...")
    adata = aggregate_batches(input_data_file_or_dir)
else:
    vprint(f"Input is a file")
    adata = sc.read_h5ad(input_data_file_or_dir)

if embedding != None:
    adata.X = adata.obsm[embedding].copy()


if not hf.is_normalized(adata):
    vprint("normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
else:
    print("Data already normalized, skipping to next step...")
vprint("computing PCA embeddings...")
sc.tl.pca(adata, svd_solver="arpack")
vprint("finding neighbors...")
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
vprint("computing leiden clusters...")
sc.tl.leiden(adata, resolution=0.1)
vprint("computing UMAP embeddings...")
sc.tl.umap(adata)

# plotting
plot_with_external_legend(adata, show=True, save=False, color="cell_type", title=f"UMAP all ductal cells for {os.path.basename(input_data_file_or_dir)}")
plot_with_external_legend(adata, show=True, save=False, color="leiden", title=f"UMAP only ductal cells colored by leiden clusters for {os.path.basename(input_data_file_or_dir)}")
"""
for annotation in annotation_columns:
    plot_with_external_legend(adata, show=False, save=True, color=annotation, title=f"UMAP colored by {annotation} for {os.path.basename(input_data_file_or_dir)}")

if os.path.isdir(input_data_file_or_dir) or ("cancer_state" in adata.obs.columns and "batch" in adata.obs.columns):
    plot_with_external_legend(adata, show=False, save=True, color="cancer_state", title=f"UMAP colored by cancer state for {os.path.basename(input_data_file_or_dir)}")
    plot_with_external_legend(adata, show=False, save=True, color="batch", title=f"UMAP colored by batch for {os.path.basename(input_data_file_or_dir)}")
"""
print(f"Output: {output_data_dir}")