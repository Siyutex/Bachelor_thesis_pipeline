# The input to this script should be a directory with preprocessed h5ad files (only filtering, no normlization, no HVG selection)
# The script aggregates all datasets of a single tissue type (cancer and non cancer) into an anndata
# Then it runs batch correction ber on that anndata using scVI, training a new model if one does not exist for the set parameters
# all of this will happen directly in gene expression space
# The results can also be plotted

import os
import sys
import scanpy as sc
import scvi # needed for batch correction
import re # needed to clean up file names
from pympler import asizeof # needed for memory profiling
import matplotlib.pyplot as plt
from global_scripts import helper_functions as hf
import json

print("batch correction script initiated")
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
    """
    Gets all cancer types present in the INPUT_DIR from file names.

    Iterates over all files in the given directory, checks if they are h5ad files and
    removes the appended numbers, "_non_cancerous" or "_cancerous" suffixes, and the "preprocessed_" prefix
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
    vprint = hf.make_vprint(verbose)

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
    list
        List of file names that have been loaded
    """
    vprint = hf.make_vprint(verbose)

    vprint("Importing batches...")
    batch_filenames = sorted([f for f in os.listdir(INPUT_DIR) if f.startswith(f"{IMPORT_PREFIX}_{cancer_type}")], key=natural_key) # sort according to numbers in suffix (e.g. 1, 2, 3)
    batch_names = [f.removesuffix(".h5ad").removeprefix(f"{IMPORT_PREFIX}_") for f in batch_filenames] # like cancer_type_cancerous_1
    cancer_state = ["non_cancerous" if "non_cancerous" in f else "cancerous" for f in batch_names] # cancerous or non_cancerous

    processed_batches = [] # list of batches that wer loaded (= used or skipped), used for metadata in executor script

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
        new_adata = sc.read_h5ad(os.path.join(INPUT_DIR, batch_file_name))
        processed_batches.append(batch_file_name)
        new_adata.obs["batch"] = batch_names[i] # assign batch label to each cell in new_adata
        new_adata.obs["cancer_state"] = state # assign cancer state to each cell in new_adata


        adata_list.append(new_adata)
        if state == "cancerous":
            cancerous_obs += new_adata.n_obs
        else:
            non_cancerous_obs += new_adata.n_obs

    # concatenate once at the end
    print("Conatenating batches...")
    if adata_list:
        adata = sc.concat(adata_list, merge="same", join="outer", index_unique="-X", fill_value=0) # join = outer so all genes are kept (not just the intersection of them)
    else:
        adata = sc.AnnData()  # empty
        print(f"No data for {cancer_type} found.")

    return adata, processed_batches

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

    print("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=max_considered_genes,      # n_top_genes is the total number of HVGs across all batches
    batch_key="batch"
    )

    n_highly_variable_genes = adata.var['highly_variable'].sum()
    print(f"Found {n_highly_variable_genes} highly variable genes across {adata.obs['batch'].nunique()} batches.")

    adata_hvg = adata[:, adata.var['highly_variable']].copy() # only keep HVGs

    #number of genes that are left
    print(f"after applying the boolean mask, there are {adata_hvg.shape[1]} genes left")
    print(f"The shape of adata_hvg is {adata_hvg.shape}")


    # only consider HVGs present in at least 30% of batches (and at very least 2 batches)
    threshold = int(0.3 * adata.obs['batch'].nunique())
    mask = adata_hvg.var.get('highly_variable_nbatches') >= threshold
    adata_hvg = adata_hvg[:, mask].copy()
    print(f"The shape of adata_hvg after thresholding is {adata_hvg.shape}")


    try:
        model = scvi.model.SCVI.load("my_scvi_model/", adata_hvg)
    except Exception as e:
        print("No scVI model found. Training a new one...")
        print(f"Model loading failed because of Exception: {e}")
        model = None
        pass

    print(model) # check if model is loaded or none

    if model is None:
        # set up model with batch information
        scvi.model.SCVI.setup_anndata(adata_hvg, batch_key="batch")

        # Train the model
        model = scvi.model.SCVI(adata_hvg)
        model.train()
        model.save("my_scvi_model/", overwrite=True)

    # Get the batch-corrected latent representation (obsm is a matrix like X where each row is a cell and each column is a feature)
    adata_hvg.obsm["X_scVI_latent"] = model.get_latent_representation()
    adata_hvg.obsm["X_scVI_corrected"] = model.get_normalized_expression()
    print(f"Shape of corrected expression matrix{adata_hvg.obsm["X_scVI_corrected"].shape}")
    print(f"Shape of latent representation{adata_hvg.obsm['X_scVI_latent'].shape}")

    return adata_hvg


def plot_with_external_legend(adata, color, **kwargs):
    """
    Plots UMAP (or other Scanpy plots) with the legend outside the figure and smaller font.
    Produces color.len() separate plots.
    
    Parameters
    ----------
    adata : AnnData
        Your annotated data object.
    color : str or list of str
        Column name(s) in adata.obs to color by.
    **kwargs : additional keyword arguments
        Any additional arguments accepted by sc.pl.umap.
    """
    # Ensure color is a list for uniform handling
    if isinstance(color, str):
        color_keys = [color]
    else:
        color_keys = color

    for key in color_keys:
        sc.pl.umap(adata, color=key, show=False, **kwargs)
        ax = plt.gca()
        # Move legend outside and shrink font
        ax.legend(
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            fontsize='small',
            title=key
        )
        plt.tight_layout()
        plt.show()

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
    corrected_adata = adata.copy()
    corrected_adata.X = corrected_adata.obsm["X_scVI_corrected"]

    print("Clustering corrected data...")
    print("Scaling...")
    sc.pp.scale(corrected_adata, max_value=10)
    print("PCA...")
    sc.pp.pca(corrected_adata, svd_solver="arpack")
    print("Neighbors...")
    sc.pp.neighbors(corrected_adata)
    print("Leiden clustering...")
    sc.tl.leiden(corrected_adata)
    print("UMAP...")
    sc.tl.umap(corrected_adata)


    # Plot to check batch correction
    print("Plotting UMAP of corrected data...")
    plot_with_external_legend(corrected_adata, color=["leiden", "batch","cancer_state"], title="UMAP of corrected data")


def plot_uncorrected_umap(adata):
    # Use corrected representation for clustering/UMAP
    """
    Plot UMAP of uncorrected data.

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
    print("Clustering uncorrected data...")

    print("Scaling...")
    sc.pp.scale(adata, max_value=10)
    print("PCA...")
    sc.pp.pca(adata, svd_solver="arpack")
    print("Neighbors...")
    sc.pp.neighbors(adata)
    print("Leiden clustering...")
    sc.tl.leiden(adata)
    print("UMAP...")
    sc.tl.umap(adata)


    # Plot to check batch correction
    print("Plotting UMAP of uncorrected data...")
    plot_with_external_legend(adata, color=["leiden", "batch","cancer_state"], title="UMAP of uncorrected data")


def get_correction_metrics():
    return 0



# execute script
cancertypes = get_cancer_types(INPUT_DIR)
for cancertype in cancertypes:
    adata, processed_batches = aggregate_batches(cancertype, max_obs_cancerous=0, max_obs_non_cancerous=50000)
    print(f"adata size in GB: {asizeof.asizeof(adata) / 1024**3}") # size of adata in GB

    adata_hvg = correct_batches(adata, 1000) 



    # save the processed data to temporary h5ad file, make relevant directory first
    final_output_path = os.path.join(OUTPUT_DIR,f"batch_corrected_HVG_{cancertype}.h5ad")
    adata_hvg.write(final_output_path)
    print("Output: " + final_output_path)

    # save metadata to temp as well
    with open(os.path.join(OUTPUT_DIR, f"metadata.json"), "w") as f:
        json.dump(processed_batches, f)
    print("Metadata: " + os.path.join(OUTPUT_DIR, f"metadata.json"))




    # plot only after saving, so I don't change the matrix anymore (eg scaling will brick future PCAs)
    # plot_corrected_umap(adata_hvg)
    # plot_uncorrected_umap(adata_hvg)


