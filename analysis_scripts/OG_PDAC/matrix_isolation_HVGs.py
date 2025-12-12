import scanpy as sc
import helper_functions as hf
import os
import warnings
import pandas as pd


def add_log1p_layer(adata):
    adata_log = sc.pp.log1p(adata, copy=True)
    new_adata = adata.copy()
    new_adata.layers["log1p"] = adata_log.X
    return new_adata


def select_genes(adata, max_considered_genes, use_log1p, genes_to_keep) -> sc.AnnData:
    
    # check whether adata.var_names are ENSG IDs or gene symbols
    if all("ENSG" in varname for varname in adata.var_names):
        varnames_format = "gene_ids"
    elif not any("ENSG" in varname for varname in adata.var_names):
        varnames_format = "gene_symbols"

    # check whether genes_to_keep are ENSG IDs or gene symbols
    if all("ENSG" in gene for gene in genes_to_keep):
        genes_to_keep_format = "gene_ids"
    elif not any("ENSG" in gene for gene in genes_to_keep):
        genes_to_keep_format = "gene_symbols"

    # remove invalid genes from genes_to_keep
    if genes_to_keep_format == "gene_ids":
        if "gene_ids" in adata.var.columns:
            invalid_genes = [g for g in genes_to_keep if g not in adata.var["gene_ids"]]
        else:
            invalid_genes = [g for g in genes_to_keep if g not in adata.var_names]
    elif genes_to_keep_format == "gene_symbols":
        if "gene_symbols" in adata.var.columns:
            invalid_genes = [g for g in genes_to_keep if g not in adata.var["gene_symbols"]]
        else:
            invalid_genes = [g for g in genes_to_keep if g not in adata.var_names]

    if invalid_genes:
        warnings.warn(f"Invalid genes found in genes_to_keep: {invalid_genes}, removing them...")
        genes_to_keep = [g for g in genes_to_keep if g not in invalid_genes]
  
    # Convert genes_to_keep if formats differ
    if genes_to_keep_format != varnames_format:
        # Case 1: converting gene_ids → gene_symbols
        if genes_to_keep_format == "gene_ids" and varnames_format == "gene_symbols":
            id_to_symbol = dict(zip(adata.var["gene_ids"], adata.var_names)) # builds a dict with the first tuple value as key and the second as value
            genes_to_keep_converted = [id_to_symbol.get(g) for g in genes_to_keep]

        # Case 2: converting gene_symbols → gene_ids
        elif genes_to_keep_format == "gene_symbols" and varnames_format == "gene_ids":
            symbol_to_id = dict(zip(adata.var["gene_symbols"], adata.var_names))
            genes_to_keep_converted = [symbol_to_id.get(g) for g in genes_to_keep]

        # Replace the list
        genes_to_keep = genes_to_keep_converted

    # set preservation mask
    preservation_mask = pd.Series(adata.var_names.isin(genes_to_keep), index=adata.var.index)

    # select HVGs, ignoring batch origin (since at this point, the data should be batch corrected)
    vprint("Selecting highly variable genes...")
    if use_log1p:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat",
            n_top_genes=max_considered_genes,      # n_top_genes is the total number of HVGs across all batches
            layer="log1p"
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=max_considered_genes,      # n_top_genes is the total number of HVGs across all batches
        )

    n_highly_variable_genes = adata.var['highly_variable'].sum()
    vprint(f"Found {n_highly_variable_genes} highly variable genes across {adata.n_obs} cells")


    final_mask = adata.var["highly_variable"] | preservation_mask
    adata_hvg = adata[:, final_mask].copy() # only keep HVGs or genes that should be preserved

    #number of genes that are left
    vprint(f"after applying the final mask (HVG | gene_to_keep), there are {adata_hvg.shape[1]} genes left")
    vprint(f"The shape of adata_hvg is {adata_hvg.shape}")

    return adata_hvg


def limit_cells(adata, isolation_dict, preservation_dict) -> sc.AnnData:
    """
    returns new adata object with only cells that fulfill all conditions in isolation_dict
    or at least one condition in preversation_dict (eg being part of a particular cnv_clade)
    """
    # isolation dict should have key = obs columns in adata, values = list of entries in that column to keep
    # in the end this function will produce the intersection of all entries (eg cell_type: ["ductal"], cancer_state:["normal","transitional"] will isolate cells that or ductal and either normal or transitional)
    
    # check that dicts correspond to valid adata.obs columns and values
    def validate_dict(dict):
        invalid_keys = set()

        for obs_column, valid_entries in dict.items():
            if obs_column in adata.obs.columns and all(entry in adata.obs[obs_column].values for entry in valid_entries): # check if column exists and needed values exist in it
                continue
            elif obs_column not in adata.obs.columns:
                warnings.warn(f"obs column {obs_column} not found in adata.obs, ignoring for isolation")
                invalid_keys.add(obs_column)
            elif not all(entry in adata.obs[obs_column].values for entry in valid_entries):
                for entry in valid_entries:
                    if entry not in adata.obs[obs_column].values:
                        warnings.warn(f"entry {entry} not found in adata.obs['{obs_column}'], ignoring column '{obs_column}' for isolation")
                        invalid_keys.add(obs_column)

        for invalid_key in invalid_keys:
            del dict[invalid_key]

    # validate isolation and preservation dicts
    validate_dict(isolation_dict)
    validate_dict(preservation_dict)

    # initialize masks
    isolation_mask = pd.Series(True, index=adata.obs.index)
    preservation_mask = pd.Series(False, index=adata.obs.index)

    # all isolation conditions must be satisfied
    for obs_column, valid_entries in isolation_dict.items():
        isolation_mask &= adata.obs[obs_column].isin(valid_entries)

    # any preservation condition may be satisfied
    for obs_column, valid_entries in preservation_dict.items():
        preservation_mask |= adata.obs[obs_column].isin(valid_entries)

    # final mask: keep cells that satisfy ALL isolation OR ANY preservation condition
    final_mask = isolation_mask | preservation_mask

    return adata[final_mask, :]


def main():

    # load adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)
    vprint(f"Adata summary:\n{adata}")
    
    print("isolating main layer...")
    adata = hf.matrix_to_anndata(adata, main_layer).copy()

    
    if add_log1p:
        print("Adding log1p layer...")
        adata = add_log1p_layer(adata).copy()
        use_log1p = True
    else:
        use_log1p = False

    # isolate cells (eg transtion state) (do before HVG to only take HVGs relevant to those cells)
    if isolation_dict != {}:
        print("limiting cells...")
        adata = limit_cells(adata, isolation_dict, preservation_dict).copy()

    # select HVGs
    if max_considered_genes != "all":
        print("selecting HVGs...")
        adata = select_genes(adata, max_considered_genes, use_log1p, genes_to_keep)
    else:
        print("Skipping HVG selection...")

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))


if __name__ == "__main__":

    input_data_file, output_data_dir, main_layer, add_log1p, max_considered_genes, genes_to_keep, isolation_dict, preservation_dict, verbose = hf.import_cmd_args(9)
    vprint = hf.make_vprint(verbose)

    main()