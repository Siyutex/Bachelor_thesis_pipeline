import scanpy as sc
import helper_functions as hf
import os
import warnings


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


def limit_cells(adata, isolation_dict) -> sc.AnnData:
    """returns new adata object with only intersection of conditions"""
    # isolation dict should have key = obs columns in adata, values = list of entries in that column to keep
    # in the end this function will produce the intersection of all entries (eg cell_type: ["ductal"], cancer_state:["normal","transitional"] will isolate cells that or ductal and either normal or transitional)
    
    vprint("copying adata")
    internal_adata = adata.copy()

    for obs_column, entry_list in isolation_dict.items():
        vprint(f"limiting {obs_column} to {entry_list}")
        if obs_column in internal_adata.obs.columns and all(entry in internal_adata.obs[obs_column].values for entry in entry_list): # check if column exists and needed values exist in it
            vprint("limiting possible")
            internal_adata = internal_adata[internal_adata.obs[obs_column].isin(entry_list)].copy()
        elif obs_column not in internal_adata.obs.columns:
            vprint("obs column not found")
            warnings.warn(f"obs column {obs_column} not found in adata.obs, skipping to next column...")
        elif not all(entry in internal_adata.obs[obs_column].values for entry in entry_list):
            for entry in entry_list:
                if entry not in internal_adata.obs[obs_column].values:
                    vprint(f"entry {entry} not found in {obs_column}")
                    warnings.warn(f"entry {entry} not found in adata.obs['{obs_column}'], skipping to next column...")

    return internal_adata


def main():

    # load adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)
    vprint(f"Adata summary:\n{adata}")
    
    print("isolating main layer...")
    adata = hf.matrix_to_anndata(adata, main_layer).copy()

    # select HVGs
    if max_considered_genes != "all":
        print("selecting HVGs...")
        select_HVGs(adata, max_considered_genes, batch_threshold)
    else:
        print("Skipping HVG selection...")

    # isolate cells (eg transtion state)
    adata = limit_cells(adata, isolation_dict)

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))


if __name__ == "__main__":

    input_data_file, output_data_dir, main_layer, max_considered_genes, batch_threshold, isolation_dict, verbose = hf.import_cmd_args(6)
    vprint = hf.make_vprint(verbose)

    main()