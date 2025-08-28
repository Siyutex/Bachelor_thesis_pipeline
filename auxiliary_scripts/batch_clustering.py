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

def aggregate_batches(directories):
    """
    Load AnnData files from given directories, treating each file as a batch.
    Returns a merged AnnData object with .obs['batch'] indicating the batch label (filename).
    """
    batches = []
    labels = []

    for d in directories:
        for fname in os.listdir(d):
            if fname.endswith(".h5ad"):   # adjust if your files are in another format
                path = os.path.join(d, fname)
                print(f"Loading {path}")
                adata = sc.read_h5ad(path)
                batches.append(adata)
                labels.append(os.path.splitext(fname)[0])  # filename w/o extension as label
    
    # compute intersection of gene names
    intersect_genes = set(batches[0].var_names)
    for ad in batches[1:]:
        intersect_genes &= set(ad.var_names)
    intersect_genes = sorted(list(intersect_genes))
    print(f"Keeping {len(intersect_genes)} genes present in all batches.")

    # subset each batch
    batches = [ad[:, intersect_genes] for ad in batches]

    # Merge all batches into one AnnData
    merged = anndata.concat(
        batches,
        join="outer",           # union of genes, fill missing with 0
        label="batch",          # new column in .obs
        keys=labels,
        fill_value=0,
        index_unique="-"
    )

    return merged


def main(directories):
    # Step 1: load and merge
    adata = aggregate_batches(directories)

    # Step 2: dimensionality reduction
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    # Step 3: plot UMAP colored by batch
    sc.pl.umap(adata, color="batch", title="UMAP colored by batch", legend_loc="on data")


if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("Usage: python batch_umap.py <dir1> <dir2> ...")
        sys.exit(1)

    dirs = sys.argv[1:]
    main(dirs)
