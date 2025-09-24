# this script will produce 2D embeddings with UMAP for visualization
# it can take any number of annotated batches as input
# if several batches are passed, it will merge them into 1 object, then cluster and plot, batch, leiden and annotations
# if one batch is passed, it will just cluster and plot, leiden and annotations 

import os
import sys
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import json
import helper_functions as hf



# import command line arguments from executor script
input_data_file_or_dir = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input file")
output_data_dir = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")
python_objects = sys.argv[3] if len(sys.argv) > 3 else print("Please provide a list with python objects. For the preprocessing subprocess, this should contain the pipeline mode as a string")

# import annotation column names
if len(python_objects) > 0:
    python_objects = python_objects.split("\n")
    for object in python_objects:
        python_objects[python_objects.index(object)] = json.loads(object.strip())
    annotation_columns = python_objects[0]
    verbose = python_objects[1]

vprint = hf.make_vprint(verbose=True)



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


def plot_with_external_legend(adata, **kwargs):
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
    plt.show()





# check if input is directory of file
if os.path.isdir(input_data_file_or_dir):
    vprint(f"Input is a directory with {len(os.listdir(input_data_file_or_dir))} files")
    vprint("aggregating batches...")
    adata = aggregate_batches(input_data_file_or_dir)
else:
    vprint(f"Input is a file")
    adata = sc.read_h5ad(input_data_file_or_dir)

vprint("normalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
vprint("computing PCA embeddings...")
sc.tl.pca(adata, svd_solver="arpack")
vprint("finding neighbors...")
sc.pp.neighbors(adata)
vprint("computing leiden clusters...")
sc.tl.leiden(adata)
vprint("computing UMAP embeddings...")
sc.tl.umap(adata)

# plotting
plot_with_external_legend(adata, color="leiden", title=f"UMAP colored by Leiden for {os.path.basename(input_data_file_or_dir)}")

for annotation in annotation_columns:
    plot_with_external_legend(adata, color=annotation, title=f"UMAP colored by {annotation} for {os.path.basename(input_data_file_or_dir)}")

if os.path.isdir(input_data_file_or_dir):
    plot_with_external_legend(adata, color="cancer_state", title=f"UMAP colored by cancer state for {os.path.basename(input_data_file_or_dir)}")
    plot_with_external_legend(adata, color="batch", title=f"UMAP colored by batch for {os.path.basename(input_data_file_or_dir)}")

