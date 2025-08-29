# this script will do PCA / UMAP and color cells according to the batch they belong to
# which will help asses how / if batch correction will be done

# 1 function to aggregate batches
# it will take directories as arguments and save each file contained therein as a batch
# merging all batches into 1 anndata object with a obs key for the batch number = filename of batch


# then just do the normal embedding

# color per batch and display as clusterplot

#!/usr/bin/env python3

import os
import sys
import scanpy as sc
import anndata

directories = [r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\colon_cancer_preprocessed",
               r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\colon_normal_preprocessed"]
sample_size = 2  # how many files to sample from each directory, to avoid memory issues

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
                batch_labels.append(fname.split("_")[1].split(".")[0])  # assuming filenames like ..._[relevant_label]. ...
                cancer_states.append("cancer" if "cancer" in d else "normal")
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


def main(directories, sample_size):
    # Step 1: load and merge
    adata = aggregate_batches(directories, sample_size)

    # Step 2: dimensionality reduction
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    # Step 3: plot UMAP colored by batch
    sc.pl.umap(adata, color="batch", title="UMAP colored by batch", legend_loc="best")


if __name__ == "__main__":

    main(directories, sample_size)
