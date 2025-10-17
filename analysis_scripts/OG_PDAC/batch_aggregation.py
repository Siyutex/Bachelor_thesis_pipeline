# aggreagate batches into one anndata object and keep all existing annotations

import os
import scanpy as sc
import re
import helper_functions as hf



# function to sort files numerically
def natural_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]


# acquire the set of cancer types present in preprocessed data from file names
def get_cancer_types(input_dir: str, input_prefix):
    """
    Gets all cancer types present in the INPUT_DIR from file names.

    Iterates over all files in the given directory, checks if they are h5ad files and
    removes the appended numbers, "_non_cancerous" or "_cancerous" suffixes, and the input_prefix prefix
    to extract the cancer type from the filename. The result is a set of unique cancer types.

    Isolated cancer types will look like: "glioblastoma", "clear_cell_renal_carcinoma"

    Parameters
    ----------
    input_dir : str
        Path to the directory containing preprocessed data as h5ad files

    Returns
    -------
    set
        Set of unique cancer types present in the preprocessed data
    """

    cancer_types = set()
    for file in os.listdir(input_dir):
        if file.endswith(".h5ad"):
            cancertype = file.removesuffix(".h5ad")
            cancertype = cancertype.removeprefix(f"{input_prefix}_")
            cancertype = re.sub(r"_\d+$", "", cancertype) # removes all numbers at the end
            if "_non_cancerous" in cancertype: cancertype = cancertype.removesuffix("_non_cancerous")
            if "_cancerous" in cancertype: cancertype = cancertype.removesuffix("_cancerous")
            cancer_types.add(cancertype) # set guarantees no duplicates
            vprint(f"Found cancer type: {cancertype}")

    vprint(f"List of unique cancer types in INPUT_DIR: {', '.join(cancer_types)}")

    return cancer_types

def aggregate_batches(cancer_type: str, input_prefix,  max_obs_cancerous: int = None, max_obs_non_cancerous: int = None):
    """
    Aggregate all batches for the given cancer type into one AnnData object. And annote the cells
    with batch and cancer state. max_obs_cancerous and max_obs_non_cancerous can be used to limit
    the total number of cells. Note that this may lead to not all batches being included.
    
    Parameters
    ----------
    cancer_type : str
        Cancer type to aggregate batches for
    max_obs_cancerous : int, optional
        Maximum number of cancerous cells to include, by default None = all batches included
    max_obs_non_cancerous : int, optional
        Maximum number of non-cancerous cells to include, by default None = all batches included

    Returns
    -------
    anndata.AnnData
        AnnData object merged based on parameters from batches for the given cancer_type
    """

    print("Importing batches...")
    batch_filenames = sorted([f for f in os.listdir(input_data_dir) if f.startswith(f"{input_prefix}_{cancer_type}")], key=natural_key) # sort according to numbers in suffix (e.g. 1, 2, 3)
    batch_names = [f.removesuffix(".h5ad").removeprefix(f"{input_prefix}_") for f in batch_filenames] # like cancer_type_cancerous_1
    cancer_state = ["non_cancerous" if "non_cancerous" in f else "cancerous" for f in batch_names] # cancerous or non_cancerous

    # make batch names shorter (cancerous needs to be first,  bcs "cancerous" is a substring of "non_cancerous")
    for i, batch_name in enumerate(batch_names):
        if "cancerous" in batch_name: batch_names[i] = "c" + batch_name.split("_")[-1]
        if "non_cancerous" in batch_name: batch_names[i] = "nc" + batch_name.split("_")[-1]


    adata_list = []
    cancerous_obs = 0
    non_cancerous_obs = 0

    vprint(f"Batch filenames: {batch_filenames}")
    vprint(f"Batch names: {', '.join(batch_names)}")
    vprint(f"Cancer state: {', '.join(cancer_state)}")

    for i, batch_file_name in enumerate(batch_filenames):
        state = cancer_state[i]

        # enforce max_obs limits
        vprint(f"Processing {state} batch {batch_file_name}")
        if state == "cancerous" and max_obs_cancerous is not None and cancerous_obs >= max_obs_cancerous:
            vprint(f"Skipping {batch_file_name} because it is {state} and {max_obs_cancerous} {state} cells have already been included.")
            continue
        if state == "non_cancerous" and max_obs_non_cancerous is not None and non_cancerous_obs >= max_obs_non_cancerous:
            vprint(f"Skipping {batch_file_name} because it is {state} and {max_obs_non_cancerous} {state} cells have already been included.")
            continue

        # read and annotate
        new_adata = sc.read_h5ad(os.path.join(input_data_dir, batch_file_name))
        new_adata.obs["batch"] = batch_names[i] # assign batch label to each cell in new_adata
        new_adata.obs["cancer_state"] = state # assign cancer state to each cell in new_adata

        adata_list.append(new_adata)
        if state == "cancerous":
            cancerous_obs += new_adata.n_obs
        else:
            non_cancerous_obs += new_adata.n_obs

    # check if var_names map to unqiue entries in each key in adata.var (if a var_name maps to more than one entry across batches, it is impossible to map the entries back to the var_name after concatenation)
    vprint("Checking if mapping var_names -> adata.var[key] is unqiue across all batches...")

    # get sets of unqiue var_names and keys across all batches
    unique_var_names = set()
    unique_keys = set()
    for adata in adata_list:
        unique_var_names.update(adata.var_names)
        unique_keys.update(adata.var.keys())

    transferability = {} # wether a key can be transfered to the concatenated adata via mapping

    # each unique_var_name has to map to exactly one entry per key across all batches (adatas in adata_list)
    for key in unique_keys:
        unique_tracker = []
        for unique_var_name in list(unique_var_names):
            key_entries = []
            for adata in adata_list: 
                if unique_var_name in adata.var_names:
                    key_entries.append(adata.var[key][unique_var_name])
            unique_key_entries = set(key_entries)
            if len(unique_key_entries) > 1:
                vprint(f"var_name {unique_var_name} has more than one unqiue entry in {key} associated with it across batches: {key_entries}")
                unique_tracker.append(False)
            else:
                unique_tracker.append(True)

        if len(unique_tracker) == len(unique_var_names) and False not in unique_tracker:
            vprint(f"All var_names across all batches correspond to exactly one entry in {key} across all batches")
            transferability[key] = True
        elif len(unique_tracker) == len(unique_var_names) and False in unique_tracker:
            vprint(f"WARNING: Not all var_names across all batches correspond to exactly one entry in {key} across all batches")
            transferability[key] =  False
        assert len(unique_tracker) == len(unique_var_names), "Logic error: tracker length mismatch"

    vprint(f"Transferability: {transferability}")

    # concatenate once at the end
    print("Concatenating batches...")
    if adata_list:
        adata = sc.concat(adata_list, merge="same", join="outer", index_unique="-X", fill_value=0) # join = outer so all genes are kept (not just the intersection of them), merge = same does not work but we remediate this afterwards by reattaching important var columns
    else:
        raise ValueError(f"No batches found for cancer type {cancer_type} in {input_data_dir}")

    # reattach adata.var[key] to aggregated adata for each key that is transferable (the concant function cannot handle it properly)
    print("Mapping var_names -> adata.var[key]...")
    for key, transferable in transferability.items():
        if transferable:
            mapping = {}
            for ad in adata_list:
                mapping.update(ad.var[key].to_dict())  # add key values pairs: "ensemb_id" : "gene_symbol"
            vprint(f"Mapping {key} to concatenated adata...")
            adata.var[key] = adata.var_names.map(mapping)
    
    vprint(f"Concatenated adata var annotations: {adata.var.keys()}")


    # test concatenated adata for duplicate var_names
    adata_dup = adata.var_names.value_counts()
    adata_dup = adata_dup[adata_dup > 1]

    if not adata_dup.empty:
        print("Duplicate genes found in concatenated AnnData (run script with verbose for details)")
        for gene, count in adata_dup.items():
            vprint(f"  {gene}: {count} entries")

    return adata


def main(input_data_dir, output_dir, max_obs_cancerous, max_obs_non_cancerous, input_prefix, output_prefix, verbose):

    cancertypes = get_cancer_types(input_data_dir, input_prefix=input_prefix) # returns a set of unique cancer types

    for cancertype in cancertypes:
        adata = aggregate_batches(cancertype, input_prefix=input_prefix, max_obs_cancerous=max_obs_cancerous, max_obs_non_cancerous=max_obs_non_cancerous)
        adata.write(os.path.join(output_dir, f"{output_prefix}_{cancertype}.h5ad"), compression="gzip")
    
    print("Output: " + output_dir)


if __name__ == "__main__":

    # import cmd args
    input_data_dir, output_dir, max_obs_cancerous, max_obs_non_cancerous, input_prefix, output_prefix, verbose = hf.import_cmd_args(4)
    vprint = hf.make_vprint(verbose)

    main(input_data_dir, output_dir, max_obs_cancerous, max_obs_non_cancerous, input_prefix, output_prefix, verbose)