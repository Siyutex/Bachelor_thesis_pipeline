# this script will produce 2D embeddings with UMAP for visualization
# it can take any number of annotated batches as input
# if several batches are passed, it will merge them into 1 object, then cluster and plot, batch, leiden and annotations
# if one batch is passed, it will just cluster and plot, leiden and annotations 

import os
import sys
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

# import command line arguments from ececutor script
input_data_file_or_dir = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input file")
output_data_dir = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")

# WIP decode from json python list object
annotations = sys.argv[3] if len(sys.argv) > 3 else print("Please provide a list with the annotation fields in adata.obs")


def aggregate_batches(directories: list[os.PathLike], sample_size: int = 1):
    """
    Load AnnData files from given directories, treating each file as a batch.
    Returns a merged AnnData object with .obs['batch'] indicating the batch label (filename).
    """
    batches = []
    batch_labels = []
    cancer_states = []

    run = 0
    for d in directories:
        print(f"Processing directory: {d}")
        

        current_dir = os.listdir(d)
        for i, fname in enumerate(current_dir):
            print(f"Index of current batch: {i}")
            if fname.endswith(".h5ad") and i+1 <= sample_size:  # adjust if your files are in another format
                path = os.path.join(d, fname)
                print(f"Loading {path}")
                adata = sc.read_h5ad(path)
                batches.append(adata)
                batch_labels.append(fname.removeprefix("preprocessed_").removesuffix(".h5ad"))  # assuming filenames like preprocessed_[relevant].h5ad
                cancer_states.append("non_cancerous" if "non_cancerous" in fname else "cancerous")
                print(batch_labels[-1])  # filename w/o extension as label
            elif not fname.endswith(".h5ad"):
                print(f"Skipping {fname}, not an .h5ad file.")
                continue
            elif i+1 >= sample_size:
                print("Sample size limit reached. Continuing to next directory.")
                break
        
    print("Cancer states of batches:", cancer_states)

    # Expand cancer states to match the number of cells in each batch
    cancer_states_expanded = []
    for state, ad in zip(cancer_states, batches):
        cancer_states_expanded.extend([state] * ad.n_obs)

    # Merge all batches into one AnnData
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

    return merged


def plot_with_external_legend(adata, **kwargs):
    # Prevent Scanpy from immediately showing the plot
    sc.pl.umap(adata, show=False, **kwargs)
    
    ax = plt.gca()
    fig = plt.gcf()
    
    # Move legend outside and shrink font
    ax.legend(
        bbox_to_anchor=(1.05, 1),  # move right of plot
        loc='upper left',
        fontsize='small',
        title=kwargs.get('color', '')  # use color key as legend title
    )
    
    # Adjust layout to avoid clipping
    plt.tight_layout()
    plt.show()



def main(directories, sample_size):
    # Step 1: load and merge
    adata = aggregate_batches(directories, sample_size)

    # Step 2: dimensionality reduction
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    # Step 3: plot UMAP colored by batch
    plot_with_external_legend(adata, color="batch", title="UMAP colored by batch")


if __name__ == "__main__":

    main(directories, sample_size)
