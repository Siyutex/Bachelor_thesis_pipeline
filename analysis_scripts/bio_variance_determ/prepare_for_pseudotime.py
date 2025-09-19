# The input to this script should be a directory with preprocessed h5ad files (only filtering, no normlization, no HVG selection)
# The script aggregates all datasets per cancer type (both cancerous and non cancerous batches) into one anndata
# it annotates batch origin and cancer state
# then log normalizes the data
# runs batch aware HVG selection
# and computes a PCA representation

import os
import sys
import scanpy as sc
import re # needed to clean up file names
import helper_functions as hf

# import command line arguments from ececutor script
INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input directory")
OUTPUT_DIR = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")

# global constants
IMPORT_PREFIX = "preprocessed" # "_" should be manually appended


# function to sort files numerically
def natural_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]


# acquire the set of cancer types present in preprocessed data from file names
def get_cancer_types(input_dir: str, verbose=False):
    vprint = hf.make_vprint(verbose)

    # get a set with all cancer type identifiers present in the filenames in the input dir
    # assumes files are coming directly from preprocessing.py
    cancer_types = set()
    for file in os.listdir(input_dir):
        if file.endswith(".h5ad"):
            cancertype = file.removesuffix(".h5ad")
            cancertype = cancertype.removeprefix(f"{IMPORT_PREFIX}_")
            cancertype = re.sub(r"_\d+$", "", cancertype) # removes all numbers at the end
            if "_non_cancerous" in cancertype: cancertype = cancertype.removesuffix("_non_cancerous")
            if "_cancerous" in cancertype: cancertype = cancertype.removesuffix("_cancerous")
            cancer_types.add(cancertype) # set guarantees no duplicates
            vprint(f"Found cancer type: {cancertype}")


    vprint(f"List of unique cancer types in INPUT_DIR: {', '.join(cancer_types)}")

    return cancer_types

def aggregate_batches(cancer_type: str, max_obs_cancerous: int = None, max_obs_non_cancerous: int = None, verbose=False):
    vprint = hf.make_vprint(verbose)

    vprint("Importing batches...")
    batch_filenames = sorted([f for f in os.listdir(INPUT_DIR) if f.startswith(f"{IMPORT_PREFIX}_{cancer_type}")], key=natural_key) # sort according to numbers in suffix (e.g. 1, 2, 3)
    batch_names = [f.removesuffix(".h5ad").removeprefix(f"{IMPORT_PREFIX}_") for f in batch_filenames] # like cancer_type_cancerous_1
    cancer_state = ["non_cancerous" if "non_cancerous" in f else "cancerous" for f in batch_names] # cancerous or non_cancerous

    # make batch names shorter (cancerous needs to be first,  bcs "cancerous" is a substring of "non_cancerous")
    for i, batch_name in enumerate(batch_names):
        if "cancerous" in batch_name: batch_names[i] = "c" + batch_name.split("_")[-1]
        if "non_cancerous" in batch_name: batch_names[i] = "nc" + batch_name.split("_")[-1]


    adata_list = []
    cancerous_obs = 0
    non_cancerous_obs = 0

    vprint(f"Batch names: {', '.join(batch_names)}")
    vprint(f"Cancer state: {', '.join(cancer_state)}")

    for i, batch_file_name in enumerate(batch_filenames):
        state = cancer_state[i]

        # enforce max_obs limits
        vprint(f"Processing batch {batch_file_name} with cancer state {state}")
        if state == "cancerous" and max_obs_cancerous is not None and cancerous_obs >= max_obs_cancerous:
            vprint(f"Skipping {batch_file_name} because it is {state} and {max_obs_cancerous} {state} cells have already been included.")
            continue
        if state == "non_cancerous" and max_obs_non_cancerous is not None and non_cancerous_obs >= max_obs_non_cancerous:
            vprint(f"Skipping {batch_file_name} because it is {state} and {max_obs_non_cancerous} {state} cells have already been included.")
            continue

        # read and annotate
        new_adata = sc.read_h5ad(os.path.join(INPUT_DIR, batch_file_name))
        new_adata.obs["batch"] = batch_names[i] # assign batch label to each cell in new_adata
        new_adata.obs["cancer_state"] = state # assign cancer state to each cell in new_adata


        adata_list.append(new_adata)
        if state == "cancerous":
            cancerous_obs += new_adata.n_obs
        else:
            non_cancerous_obs += new_adata.n_obs

    # concatenate once at the end
    print(f"Concatenating batches...")
    if adata_list:
        adata = sc.concat(adata_list, merge="same", join="outer", index_unique="-X", fill_value=0) # join = outer so all genes are kept (not just the intersection of them)
    else:
        adata = sc.AnnData()  # empty
        print(f"No data for {cancer_type} found.")

    return adata


def prepare_adata(adata, max_considered_genes: int = None, reduced_dim: int = 50, verbose=False): 
    vprint = hf.make_vprint(verbose=True)

    # log normalize adata
    vprint("Log normalizing adata...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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
    threshold = int(0.3 * adata.obs['batch'].nunique())
    mask = adata_hvg.var.get('highly_variable_nbatches') >= threshold
    adata_hvg = adata_hvg[:, mask].copy()
    vprint(f"The shape of adata_hvg after thresholding is {adata_hvg.shape}")

    # run pca on hvg isoltated adata
    sc.pp.pca(adata_hvg, n_comps=reduced_dim, svd_solver="arpack")

    return adata_hvg









# ------ main part ------
print(f"Getting list of cancer types in {INPUT_DIR}")
cancer_types = get_cancer_types(INPUT_DIR, verbose=True)

for cancer_type in list(cancer_types):
    print(f"Aggreationg batchtes for {cancer_type}")
    adata = aggregate_batches(cancer_type, max_obs_cancerous=50000, max_obs_non_cancerous=50000, verbose=True)
    print(f"Preparing adata for {cancer_type}")
    adata = prepare_adata(adata, max_considered_genes=3000, verbose=True)
    final_output_path = os.path.join(OUTPUT_DIR,f"prepared_{cancer_type}.h5ad")
    adata.write(final_output_path, compression="gzip", )
    print(f"Adata for {cancer_type} saved to {final_output_path}")
    
# since we have multiple, files we just forward the entire directory to the executor script
# handle this in executor with adjusted behavior
print("Output: " + OUTPUT_DIR)