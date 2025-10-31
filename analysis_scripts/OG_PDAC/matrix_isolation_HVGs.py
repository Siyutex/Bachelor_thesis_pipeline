import scanpy as sc
import numpy as np
import helper_functions as hf
import os


def select_HVGs(adata, max_considered_genes, batch_threshold):
    
    # do batch aware HVG selection
    vprint("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=max_considered_genes,      # n_top_genes is the total number of HVGs across all batches
        batch_key="batch"
    )

    n_highly_variable_genes = adata.var['highly_variable'].sum()
    vprint(f"Found {n_highly_variable_genes} highly variable genes across {adata.obs['batch'].nunique()} batches.")

    adata_hvg = adata[:, adata.var['highly_variable']].copy() # only keep HVGs

    #number of genes that are left
    vprint(f"after applying the boolean mask, there are {adata_hvg.shape[1]} genes left")
    vprint(f"The shape of adata_hvg is {adata_hvg.shape}")


    # only consider HVGs present in at least 30% of batches (and at very least 2 batches)
    threshold = int(batch_threshold * adata.obs['batch'].nunique())
    mask = adata_hvg.var.get('highly_variable_nbatches') >= threshold
    adata_hvg = adata_hvg[:, mask].copy()
    vprint(f"The shape of adata_hvg after thresholding is {adata_hvg.shape}")

    adata = adata_hvg.copy()


def main():

    # load adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)
    vprint(f"Adata summary:\n{adata}")
    
    adata = hf.matrix_to_anndata(adata, main_layer).copy()

    # select HVGs
    if max_considered_genes != "all":
        print("selecting HVGs...")
        select_HVGs(adata, max_considered_genes, batch_threshold)
    else:
        print("Skipping HVG selection...")

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))


if __name__ == "__main__":

    input_data_file, output_data_dir, main_layer, max_considered_genes, batch_threshold, verbose = hf.import_cmd_args(5)
    vprint = hf.make_vprint(verbose)

    main()