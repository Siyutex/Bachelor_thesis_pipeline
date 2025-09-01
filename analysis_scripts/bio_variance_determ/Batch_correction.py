# this scirpt is meant to be run as a subprocess of Biological_variance_pipeline_executor.py
# it will aggregate all datasets of a single tissue type (cancer and non cancer) into an anndata
# then annote all cells with lineage based on marker genes
# then run batch correction ber lineage (this will hopefully reduce batch effects aslo between cancer and non cancer cells)
# all of this will happen directly in gene expression space

import os
import sys
import tempfile
import scanpy as sc

# import command line arguments from ececutor script
input_data_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input directory")
output_data_path = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")


# acquire the set of cancer types present in preprocessed data from file names
def get_cancer_types(input_data_path: str) -> set:
    """
    Gets all cancer types present in the preprocessed data from the file names.

    Iterates over all files in the given directory, checks if they are h5ad files and
    removes the "_non_cancerous" or "_cancerous" suffixes and the "preprocessed_" prefix
    to extract the cancer type from the filename. The result is a set of unique cancer types.

    Parameters
    ----------
    input_data_path : str
        Path to the directory containing the preprocessed data as h5ad files

    Returns
    -------
    set
        Set of unique cancer types present in the preprocessed data
    """
    cancer_types = set()
    for file in os.listdir(input_data_path):
        if file.endswith(".h5ad"):
            file.removesuffix(".h5ad")
            file.removeprefix("preprocessed_")
            if "non_canerous" in file: file.removesuffix("_non_canerous")
            if "cancerous" in file: file.removesuffix("_cancerous")
            cancer_types.add(file) # set guarantees no duplicates
    return cancer_types

def aggregate_batches(cancer_type: str, max_obs_cancerous: int = None, max_obs_non_cancerous: int = None):
    # get list of all h5ad files for the current cancer type + annotations
    batch_filenames = [f for f in os.listdir(input_data_path) if f.startswith(f"preprocessed_{cancer_type}")]
    batch_names = [f.removesuffix(".h5ad").removeprefix("preprocessed_") for f in batch_filenames]
    cancer_state = ["cancerous" if "cancerous" in f else "non_cancerous" for f in batch_names]

    adata_list = []
    cancerous_obs = 0
    non_cancerous_obs = 0

    for i, batch_file_name in enumerate(batch_filenames):
        state = cancer_state[i]

        # enforce max_obs limits
        if state == "cancerous" and max_obs_cancerous is not None and cancerous_obs >= max_obs_cancerous:
            continue
        if state == "non_cancerous" and max_obs_non_cancerous is not None and non_cancerous_obs >= max_obs_non_cancerous:
            continue

        # read and annotate
        new_adata = sc.read_h5ad(os.path.join(input_data_path, batch_file_name))
        new_adata.obs["batch"] = batch_names[i]
        new_adata.obs["cancer_state"] = state

        adata_list.append(new_adata)
        if state == "cancerous":
            cancerous_obs += new_adata.n_obs
        else:
            non_cancerous_obs += new_adata.n_obs

    # concatenate once at the end
    if adata_list:
        adata = sc.concat(adata_list, merge="same")
    else:
        adata = sc.AnnData()  # empty
        print(f"No data for {cancer_type} found.")

    return adata



# execute script
cancertypes = get_cancer_types(input_data_path)
for cancertype in cancertypes:
    aggregate_batches(cancertype)