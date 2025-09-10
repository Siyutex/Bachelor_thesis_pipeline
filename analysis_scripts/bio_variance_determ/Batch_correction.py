# this scirpt is meant to be run as a subprocess of Biological_variance_pipeline_executor.py
# it will aggregate all datasets of a single tissue type (cancer and non cancer) into an anndata
# then annote all cells with lineage based on marker genes
# then run batch correction ber lineage (this will hopefully reduce batch effects aslo between cancer and non cancer cells)
# all of this will happen directly in gene expression space

import os
import sys
import tempfile
import scanpy as sc
import scvi # needed for batch correction
import re # needed to clean up file names
from pympler import asizeof # needed for memory profiling

print("batch correction script initiated")
# import command line arguments from ececutor script
input_data_dir = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the input directory")
output_data_dir = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to the output directory")



# function to sort files numerically
def natural_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]


# acquire the set of cancer types present in preprocessed data from file names
def get_cancer_types(input_data_dir: str) -> set:
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
    for file in os.listdir(input_data_dir):
        if file.endswith(".h5ad"):
            cancertype = file.removesuffix(".h5ad")
            cancertype = cancertype.removeprefix("preprocessed_")
            cancertype = re.sub(r"_\d+$", "", cancertype) # removes all numbers at the end
            if "non_cancerous" in cancertype: cancertype = cancertype.removesuffix("_non_cancerous")
            if "cancerous" in cancertype: cancertype = cancertype.removesuffix("_cancerous")
            cancer_types.add(cancertype) # set guarantees no duplicates

    print(f"Found cancer types: {', '.join(cancer_types)}")

    return cancer_types

def aggregate_batches(cancer_type: str, max_obs_cancerous: int = None, max_obs_non_cancerous: int = None):
    # get list of all h5ad files for the current cancer type + annotations
    """
    Aggregate all batches for the given cancer type into one AnnData object.
    And annote the cells with batch and cancer state.
    
    Parameters
    ----------
    cancer_type : str
        Cancer type to aggregate batches for
    max_obs_cancerous : int, optional
        Maximum number of cancerous cells to include, by default None
    max_obs_non_cancerous : int, optional
        Maximum number of non-cancerous cells to include, by default None

    Returns
    -------
    anndata.AnnData
        AnnData object containing all batches for the given cancer type
    """
    batch_filenames = sorted([f for f in os.listdir(input_data_dir) if f.startswith(f"preprocessed_{cancer_type}")], key=natural_key)
    batch_names = [f.removesuffix(".h5ad").removeprefix("preprocessed_") for f in batch_filenames] # like cancer_type_1
    cancer_state = ["non_cancerous" if "non_cancerous" in f else "cancerous" for f in batch_names] # cancerous or non_cancerous

    adata_list = []
    cancerous_obs = 0
    non_cancerous_obs = 0

    print(f"Batch filenames: {batch_filenames}")
    print(f"Batch names: {', '.join(batch_names)}")
    print(f"Cancer state: {', '.join(cancer_state)}")

    for i, batch_file_name in enumerate(batch_filenames):
        state = cancer_state[i]

        # enforce max_obs limits
        print(f"Processing {state} batch {batch_file_name}")
        if state == "cancerous" and max_obs_cancerous is not None and cancerous_obs >= max_obs_cancerous:
            print(f"Skipping {batch_file_name} because it is {state} and {max_obs_cancerous} {state} cells have already been included.")
            continue
        if state == "non_cancerous" and max_obs_non_cancerous is not None and non_cancerous_obs >= max_obs_non_cancerous:
            print(f"Skipping {batch_file_name} because it is {state} and {max_obs_non_cancerous} {state} cells have already been included.")
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

    # concatenate once at the end
    if adata_list:
        adata = sc.concat(adata_list, merge="same", join="outer", index_unique="-X", fill_value=0) # join = outer so all genes are kept (not just the intersection of them)
    else:
        adata = sc.AnnData()  # empty
        print(f"No data for {cancer_type} found.")

    return adata

def correct_batches(adata, max_considered_genes: int = None):
    """
    Correct batches using scVI. Adds a batch-corrected latent representation to `adata.obsm['X_scVI']`.
    Does internal normalization, so feed raw (filtered) counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data with batch information in `adata.obs['batch']`.

    Returns
    -------
    adata : anndata.AnnData
        Input data with batch-corrected latent representation in `adata.obsm['X_scVI']`.
    """

    # do batch aware HVG selection
    B = adata.obs['batch'].nunique() # number of batches
    if max_considered_genes is not None: HVGs_per_batch = max(100, max_considered_genes // B)
    else: HVGs_per_batch = 100

    sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=HVGs_per_batch,      # per-batch top
    batch_key="batch"
    )

    adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # only consider HVGs present in at least 30% of batches (and at very least 2 batches)
    threshold = max(2, int(0.3 * B))
    mask = adata.var.get('highly_variable_nbatches' / adata.obs['batch'].nunique) >= threshold
    adata_hvg = adata[:, mask].copy()

    # set up model with batch information
    scvi.model.SCVI.setup_anndata(adata_hvg, batch_key="batch")

    # Train the model
    model = scvi.model.SCVI(adata_hvg)
    model.train()

    # Get the batch-corrected latent representation
    adata_hvg.obsm["X_scVI_latent"] = model.get_latent_representation()

    return adata_hvg


def plot_corrected_umap(adata):
    # Use corrected representation for clustering/UMAP
    """
    Plot UMAP after batch correction.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data with batch information in `adata.obs['batch']` and
        batch-corrected latent representation in `adata.obsm['X_scVI']`.

    Returns
    -------
    None

    Notes
    -----
    This function creates a UMAP plot of the data in the corrected latent
    space. The plot is colored by batch and cell type.
    """
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    # Plot to check batch correction
    sc.pl.umap(adata, color="batch", legend_loc="best")




# execute script
cancertypes = get_cancer_types(input_data_dir)
for cancertype in cancertypes:
    adata = aggregate_batches(cancertype, max_obs_cancerous=10000, max_obs_non_cancerous=10000)
    print(f"adata size in GB: {asizeof.asizeof(adata) / 1024**3}") # size of adata in GB

    adata_hvg = correct_batches(adata) 
    plot_corrected_umap(adata_hvg)


    # save the processed data to temporary h5ad file, make relevant directory first
    final_output_path = os.path.join(output_data_dir,f"batch_corrected_{cancertype}.h5ad")
    adata.write(final_output_path)
    print("Output: " + final_output_path)