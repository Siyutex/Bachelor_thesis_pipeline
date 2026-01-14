# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients
# change RAW_DATA_DIR and outcome_storage for each run to reflect new input and what outcomes should be permanenlty saved

import os # needed for file and directory operations
from enum import Enum
import shutil # needed for file storage operations
import tempfile # needed for temporary file operations
import sys # needed to exit the program 
from dataclasses import dataclass
import helper_functions as hf
from typing import Literal
import numpy as np



#check if the required directories exist, if not create them
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data")):  # check if data directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data"))  # create data directory if it does not exist
    print("Created Data directory. Please add subdirectories for each sample containing the mtx and tsv outputs of the 10x genomics pipeline")  # inform user to add subdirectories with required files
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "output_storage")):  # check if storage directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "output_storage"))  # create storage directory if it does not exist
    print("Created output_storage directory. Intermediate output files will be saved here, if specified in OUTCOME_STORAGE.")  # inform user that intermediate outputs can be saved here
if not os.path.exists(os.path.join(tempfile.gettempdir(),"python")):
    os.makedirs(os.path.join(tempfile.gettempdir(),"python"))
    print("created \"python\" directory with tempfile library to store temporary pipeline files. These files will be deleted when purgetempfiles() is called.")



# path constants
SCRIPT_DIR = os.path.dirname(__file__)  # directory where this script is located
# list of directories (see choose_pipeline_mode for valid structures for each entry)
RAW_DATA_DIRS = [    
                os.path.join(SCRIPT_DIR, "..", "..", "Data","OG_data", "manual_cell_IDs_NCBI","PDAC_cancerous"),
                os.path.join(SCRIPT_DIR, "..", "..", "Data","OG_data","manual_cell_IDs_NCBI","PDAC_non_cancerous")
                ]
OUTPUT_STORAGE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data", "output_storage")  # directory for optional permanent storage of indermediate subprocess outputs
TEMP_DIR = os.path.join(tempfile.gettempdir(),"python") # directory for storage of temporary pipeline files
AUX_DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "auxiliary_data")


# variable to determine what intermediate files should be saved permanently, one key per script
OUTCOME_STORAGE = {
    "Preprocessing.py": True,
    "Cell_type_annotation.py": True,
    "Clustering.py": False,
    "Batch_correction.py": True,
    "GRN_edge_inference.py": False,
    "GRN_rule_inference.py": False,
    "GRN_simulation.py": False,
    "infer_CNV.py": True,
    "pseudotime_inference.py": True,

    "prepare_for_pseudotime.py": False,
    "Variance.py": False
}

# classes
class pipeline_mode(Enum):
    MTX_TSVs_in_subfolders = 1  #10x genomics, format
    compressed_MTX_TSVs_in_subfolders = 2 #GDC format
    dot_matrix_files = 3 #cancerSCEM format
    h5ad_files = 4,
    tenx_h5_files = 5
    NO_MODE_CHOSEN = None

# functions
def choose_pipeline_mode(raw_data_dir):
    """Inspects the structure of the provided RAW_DATA_DIR and returns the appropraote
    pipeline_mode enum value."""
    
    mode = pipeline_mode.NO_MODE_CHOSEN

    # check if RAW_DATA_DIR has subdirectories (to avoid errors), if so check their structure
    if any(os.path.isdir(os.path.join(raw_data_dir, folder)) for folder in os.listdir(raw_data_dir)):

        if any(any(file.endswith(".mtx") for file in os.listdir(os.path.join(raw_data_dir,folder))) for folder in os.listdir(raw_data_dir)):
            # if there is an mtx in any subfolder of RAW_DATA_DIR, we assume that it is full of folders with mtx + 2 tsvs (10x genomics format)
            mode = pipeline_mode.MTX_TSVs_in_subfolders

        if any(any(file.endswith("matrix.tar.gz") for file in os.listdir(os.path.join(raw_data_dir,folder))) for folder in os.listdir(raw_data_dir)):
            # if there is such a structure, but with .gz compression (GDC format)
            mode = pipeline_mode.compressed_MTX_TSVs_in_subfolders
    elif any(file.endswith("counts.matrix.tsv.gz") for file in os.listdir(raw_data_dir)):
        # if the RAW_DATA_DIR directly contains .matrix files (cancerSCEM format)
        mode = pipeline_mode.dot_matrix_files
    elif any(file.endswith(".h5ad") for file in os.listdir(raw_data_dir)):
        mode = pipeline_mode.h5ad_files
    elif any(file.endswith(".h5") for file in os.listdir(raw_data_dir)):
        mode = pipeline_mode.tenx_h5_files

    if mode == pipeline_mode.NO_MODE_CHOSEN:
        raise ValueError("Could not determine pipeline mode. Please check the structure of the provided RAW_DATA_DIR. It should either contain subdirectories with mtx and tsv files (10x genomics format), subdirectories with compressed mtx and tsv files (GDC format) or directly .matrix files (cancerSCEM format).")

    print(f"Pipeline mode: {mode.name}")
    return mode


def purge_tempfiles(a=None, b=None): # a and b are needed for signal handler, but not used
    """
    Deletes all files and folders in the system temp directory, which is defined by
    the TEMP_DIR variable. This function is intended to be used as a signal handler
    to clean up temporary files when the program is interrupted.

    Parameters:
    a (None): not used
    b (None): not used
    
    Returns:
    None
    """
    for filename in os.listdir(TEMP_DIR):
        file_path = os.path.join(TEMP_DIR, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)       # delete file or symlink
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)   # delete folder and contents
        except Exception as e:
            print(f"Failed to delete {file_path}: {e}")

        
def request_stop(signum, frame):
    """
    Signal handler to request stopping the main loop. Also purges tempfiles.

    Parameters:
    signum (int): signal number
    frame (frame): current stack frame

    Returns:
    None
    """
    
    print(f"Caught signal {signum}, will stop soon pipeline execution...")
    purge_tempfiles()
    print("Temporary files purged.")
    sys.exit(0)


def reduce_data(
        input_data_file: str, 
        input_prefix: str,
        layers_to_remove: list[str],
        save_output: bool = False,
        output_prefix: str = "reduced",
        verbose: bool = False) -> list[str]:
    """ 
    Reduce the size of an h5ad file by removing unwanted layers / obsm matrices or adata.X.

    Input should be an h5ad file.

    Requires the following annotations to be present:
        - layers / obsm that are passed to layers_to_remove

    Outputs a gzip compressed h5ad file with reduced layers / obsm matrices
    Output files are named {output_prefix}_{basename}.h5ad.

    Adds no annotations to adata.

    Parameters:
        input_data_file (str): path to h5ad file.
        input_prefix (str): prefix of input file names, must match or will cause error.
        layers_to_remove (list[str]): list of layers / obsms to remove. Pass "X" to remove adata.X.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/reduced. Defaults to False.
        output_prefix (str, optional): prefix for output file names. Defaults to "reduced".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "reduced"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "reduced"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "reduced")
    output_temp_dir = os.path.join(TEMP_DIR, "reduced")

    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Reducing {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Data_reduction.py"), input_data_file, output_temp_dir, [layers_to_remove, verbose])

    # rename output file
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def subsample_cells(file_path, fraction, n_samples, output_dir):
    import scanpy as sc
    """
    Produce new adata with layer as X, then subsample to fraction of cells without replacement
    """

    # make sure output_dir exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # create isolated adata
    adata = sc.read(file_path)

    # create boolean mask
    n_true = int(fraction * adata.shape[0]) # number of cells to keep
    mask = np.zeros(adata.shape[0], dtype=bool)
    mask[:n_true] = True # array with size n_obs, but ordered
    
    # subsample
    for i in range(n_samples):
        print(f"Subsampling {i}")
        np.random.shuffle(mask) # shuffle to keep a random set of cells
        ss_adata = adata[mask,:].copy()
        ss_adata.write(os.path.join(output_dir, f"sample_{i}.h5ad"), compression="gzip")
        del ss_adata # free up memory


def cluster_and_plot(
        # necessary arguments
        modules: list[Literal["projections+DEGs","projections", "pseudotime_vs_cnv", "phylogenetic_tree", "gene_expression_vs_obs"]],

        # general arguments
        input_data_file: str = "", 
        selection_criteria: dict[str,list[str]] = {},
        obs_annotations: list[str] = None,

        # UMAP / PCA plotting arguments
        layers: list[str] = ["X"],
        projection: Literal["UMAP", "PCA"]  = "UMAP",
        marker_file_path: str = None,
        root_cell_idx: int = None,

        # Phylogenetic tree arguments
        tree_file: str = None,
        target_circumference: float =  0.9,
        sort_order: Literal["lowest_width_first", "lowest_depth_first"] = "lowest_width_first",

        # genes vs obs arguments
        x_axis: str = None,
        genes: list[str] = None,

        # auxiliary arguments
        show: bool = True,
        save_output: bool = False,
        output_storage_subdir = "",
        verbose: bool = False) -> None:
    
    """
    Runs plotting modules passed in modules list.

    Current modules include:
        - projections+DEGs, which will produce UMAP or PCA plots computed from the passed layer
          and colored by the passed obs annotations + leiden clusters. Also computes DEGs for categorical obs annotations + leiden clusters.
          This module uses the following arguments: input_data_file, selection_criteria, obs_annotations, layers, projection, marker_file_path, root_cell_idx
        - projections, identical to projections+DEGs, but does not compute DEGs
        - pseudotime_vs_cnv, which will produce a plot scatterplot of cells, with cnvs on the y axis and pseudotime on the x axis.
          This module uses the following arguments: input_data_file
        - phylogenetic_tree, which will produce a phylogenetic tree for the passed tree_file. 
          Can use obs annotations from anndata to draw colored annoation rings around the tree.
          This module uses the following arguments: tree_file, target_circumference, sort_order, (input_data_file and obs_annotations for colored annotation rings)
        - gene_expression_vs_obs, whill will produce a scatter plot, with the expression level of a list of genes on the y axis and a numerical obs annotation as the x axis (sorted low to high)
          
    Input should be an h5ad file with obs annotations for each cell that should be used for clustering.

    Annotations that need to be present:
        - a column adata.obs for each key in selection_criteria (and the columns need to contain at least one instance of each string in the list that is the corresponding value)
        - a column in adata.obs for each passed obs annotation
        - the specified layers need to exist in adata.layers or adata.obsm
        - if modules includes "pseudotime_vs_cnv", adata.obs["dpt_pseudotime"] (from pseudotime_inference.py), adata.obs["cnv_score"], adata.obs["summed_cnvs] (from infer_CNV.py)
        - if x_axis is not None, adata.obs[x_axis]
        - if genes is not None, adata.var_names must contain all entries in genes
        
    Outputs png files of plots and text files with DEGs if save_output is True. Shows plots if show is True.

    Does not add annotations, as it does not return an h5ad file.

    Parameters:
        modules (list[Literal["projections", "pseudotime_vs_cnv"]]): list of modules in plotting scirpt to run.
        
        input_data_file (str): path to h5ad file with obs annotations for each cell.
        selection_criteria (dict[str,list[str]]): dict of obs annotations and list of entries in that obs annotation to limit cells to (cells have to fulfill all criteria to be kept)
        obs_annotations (list[str], optional): list of obs annotations to cluster. eg ["cell_type", "pseudotime", "cnv_score", ...].
        
        layers (str, optional): Anndata layer or obsm keys to use as base for UMAP / PCA projection. Eg. ["X_cnv", "X_scvi_corrected", ...]. Defaults to ["X"], meaning adata.X will be used. 
        projection (literal ["UMAP", "PCA"], optional): projection to use for clustering. Defaults to "UMAP".
        marker_file_path (str, optional): path to file with marker genes for each cell type. DEG analysis will not yield matrix plots if not provided.
        root_cell_idx (int, optional): Index of the root cell of pseudotime. Will be highlighted in pseudotime colored plot.
        
        tree_file (str, optional): path to newick tree file. Required to plot a phylogenetic tree.
        target_circumference (float, optional): Target circumference of the phylogenetic tree. 0.9 (90% of a circle = 2*pi*0.9)
        sort_order (Literal["lowest_width_first", "lowest_depth_first"], optional): Sort order of leaves in phylogenetic tree. Defaults to "lowest_width_first". 
            lowest_width_first: Clades with the lowest number of terminal nodes are visited first ("narrowest" clades first).
            lowest_depth_first: Clades with the shortest path to a preterminal node (a node that only has terminal nodes as children) are visited first ("shallow" clades first).
            (note that the "lowest_depth_first" option may also choose to visit a generally deep clade, which happens to have a single shallow child clade, first)

        x_axis (str, optional): obs annotation to use for x axis of scatterplot. Defaults to None.
        genes (list[str], optional): list of genes to use for y axis of scatterplot. Defaults to None.    
        
        show (bool, optional): whether to show the plots. Defaults to True.
        save_output (bool, optional): whether to save the plots permanently. Defaults to False.
        output_storage_subdir (str, optional): subdirectory of "plots" to store this run's plots in.
        verbose (bool, optional): whether to print verbose output. Defaults to False.

    Returns:
        None, but saves files with plots if save_output is True
    """
    # sanity checks
    if show == False and save_output == False:
        raise ValueError("At least one of show or save_output must be True")
    for value in selection_criteria.values(): # make sure values are actual lists and not strings (lists of characters)
        assert isinstance(value, list), f"selection_criteria value {value} is not a list"
    assert isinstance(modules, list), f"modules is not a list"
    for module in modules:
        assert module in ["projections", "projections+DEGs", "pseudotime_vs_cnv", "phylogenetic_tree", "gene_expression_vs_obs"], f"module {module} is not valid"

    # adjust parameters
    target_circumference = target_circumference * 2 * np.pi

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have relevant folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "plots"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "plots"), exist_ok=True)

    if output_storage_subdir != "":
        os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "plots", output_storage_subdir), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "plots", output_storage_subdir)
    output_temp_dir = os.path.join(TEMP_DIR, "plots")

    # run script and assign path to temporary output file
    print(f"Clustering and plotting: {input_data_file}")
    hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Plotting.py"), input_data_file, output_temp_dir, [modules, selection_criteria, obs_annotations, layers, projection, marker_file_path, root_cell_idx, tree_file, target_circumference, sort_order, x_axis, genes, show, verbose, save_output])

    # naming happens in subprocess (relies on knowing which obs column was used)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        for file in os.listdir(output_temp_dir):
            shutil.copy(os.path.join(output_temp_dir, file), os.path.join(output_storage_dir, file))

    # outputs will not be used programatically so no need to return list with output file paths
    # this is why we also don't use input and output prefixes
    return None


@dataclass
class FilteringParameters:
    """
    Configuration for filtering low-quality cells and genes from single-cell data.
    If any mitochondrial filtering parameter is None, the corresponding filtering step will not be carried out. 

    Attributes:
        min_n_genes_percentile: 
            Cells with fewer expressed genes than this percentile are removed.
        min_n_cells_percentage (float): 
            Genes expressed in fewer than this fraction of cells are removed.
        min_n_UMIs_percentile:
            Cells with fewer UMI counts than this percentile are removed.
        max_n_MADs:
            Cells with mitochondrial percentages greater than (median + this many MADs) are removed.
        max_mito_percentage (float):
            Cells with mitochondrial percentages greater than this are removed.
        expected_doublet_percentage:
            Fraction of cells expected to be doublets, used by Scrublet. 
            See the sequencing device manufacturer’s recommendations.
    """

    # most default values chosen based on manual inspection of plots of OG_PDAC data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966#:~:text=Summary%20Pancreatic%20ductal%20adenocarcinoma%20,plot%20to%20predict%20the%20overall
    min_n_genes_percentile: int = 10
    min_n_cells_percentage: float = 0.01
    min_n_UMIs_percentile: int = 10
    max_n_MADs: int = None
    max_mito_percentage: float = 0.10 # default value chosen according to https://pmc.ncbi.nlm.nih.gov/articles/PMC8599307/?utm_source=chatgpt.com
    expected_doublet_percentage: float = 0.024

def preprocess_data(
        raw_data_dir: str,
        pipeline_mode: pipeline_mode, 
        filtering_params: FilteringParameters = FilteringParameters(),
        use_ensembl_ids: bool = True, 
        save_output: bool = False,
        output_prefix: str = "preprocessed", 
        verbose: bool = False,
        ) -> list[str]:
    """
    Loops through raw_data_dir and runs Preprocessing.py on each file / folder contained in it.
    
    Input should be a directory with containing scRNAseq samples in any supported format (10x genomics, GDC, cancerSCEM).

    Outputs a gzip compressed h5ad file for each sample, filtered for low-quality cells and genes and with doublets removed.
    Output files are named preprocessed_{raw_data_dir}_{index}.h5ad.

    Annotations added to adata.obs: [numpy.int64:"n_genes", numpy.float32:"n_counts", numpy.float32:"pct_counts_mito", numpy.float64:"doublet_score", numpy.bool:"predicted_doublet"] (number of genes in a cell,
    number of UMI counts in a cell, mitochondrial gene percentage in a cell, doublet score (used by Scrublet to determine whether a cell is predicted to be a doublet), whether a cell is predicted to be a doublet

    Annotations added to adata.var: [numpy.int64:"n_cells", numpy.bool:"mito"] (number of cells in which a gene is expressed, whether a gene is mitochondrial)

    Conditional annotations added to adata.var: [str:"gene_symbols" or str:"gene_ids"] (gene symbols or ensembl ids, if use_ensembl_ids is True, ensembl ids will be adata.var_names and gene symbols will be adata.var["gene_symbols"] and vice versa)

    Parameters:
        raw_data_dir (str): path to raw data directory. This should be a directory that contains datasets in any supported format.
            Currently supported formats are: 10x genomics, GDC, cancerSCEM. Example: .../glioma/dataset1, dataset2 ...
        pipeline_mode (pipeline_mode): pipeline mode enum value, run choose_pipeline_mode() to get this value.
            Used to determine how to load data for different formats from raw_data_dir.
        filtering_params (FilteringParameters, optional): configuration for filtering low-quality cells and genes from single-cell data.
            Defaults to FilteringParameters() default values. See class definition for details.
        use_ensebml_ids (bool, optional): whether to use ensembl ids for gene names. Defaults to True.
            Not all input datasets have ensembl ids (eg cancerSCEM). If false, gene symbols will be used,
            which may lead to duplicate gene names and issues downstream.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/preprocessed.
            Defaults to False.
        output_prefix (str, optional): prefix for output file names. Defaults to "preprocessed".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        output_file_list (list[str]): list of file paths to output files
    """

    filtering_params_list = [
        filtering_params.min_n_genes_percentile,
        filtering_params.min_n_cells_percentage,
        filtering_params.min_n_UMIs_percentile,
        filtering_params.max_n_MADs,
        filtering_params.max_mito_percentage,
        filtering_params.expected_doublet_percentage,
    ]

    # check if OUTCOME_STORAGE_DIR and TEMP_DIR have preprocessed folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "preprocessed"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "preprocessed")
    output_temp_dir = os.path.join(TEMP_DIR, "preprocessed")

    # assign output file list
    output_file_list = []

    # assign datatype
    datatype = str(pipeline_mode.name)

    i = 0 # iterator for file naming
    for element in os.listdir(raw_data_dir): # element can be a file or folder, data loading handled by Preprocessing.py
        print("Currently preprocessing: " +  element)
        temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, element), output_temp_dir, [datatype, filtering_params_list, use_ensembl_ids, verbose])

        # rename output file
        os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
        temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(raw_data_dir)}_{i}.h5ad")

        # add output file to output_file_list
        output_file_list.append(temp_output_path)

        # if specified, permanently store a copy of the temporary output file
        if save_output == True:
            print(f"Saving {temp_output_path} to {output_storage_dir}")
            shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

        i += 1

    return output_file_list


def annotate_cell_types(
        input_data_dir: str, 
        use_ensembl_ids: bool, 
        marker_file_path: str, 
        negative_marker_file_path: str = None,
        model: Literal["z_score","cellassign"] = "cellassign",
        cutoff_unsure: float = 0.8, 
        cutoff_other: float = -0.2, 
        save_output: bool = False, 
        input_prefix: str = "preprocessed",
        output_prefix: str = "cell_type_annotated",
        verbose: bool = False) -> list[str]:
    """ 
    Loops through input_data_dir and runs Cell_type_annotation.py on each contained h5ad file.
    Depending on the model paramter, either z_score or cellassign will be used.
    
    Z_score gets the z score of all genes after normalization, then sums up the z scores for
    set of marker corresponding to a cell type, cutoff unsure and cutoff other are then used to assign identities.
    
    Cellassign uses cellassign from scvitools to predict cell types: https://docs.scvi-tools.org/en/1.3.3/user_guide/models/cellassign.html
    A model is train for each intitial run with a unique file basename, and then saved.
    The model is reused for all subsequent runs with the same file basename to guarantee reproducibility.
    This assumes the same file basename contains the exact same data.

    e.g.:
    input_data_dir = "/path/to/preprocessed" -> contains files "PDAC_cancerous_1.h5ad", "PDAC_cancerous_2.h5ad", "PDAC_non_cancerous_1.h5ad", "PDAC_non_cancerous_2.h5ad"
    basenames will be "PDAC_cancerous_1.h5ad", "PDAC_cancerous_2.h5ad", "PDAC_non_cancerous_1.h5ad", "PDAC_non_cancerous_2.h5ad"
    one model is saved for each basename (can only be overriden if model file is manually deleted)
    reproducibility not guaranteed, if file contents change between runs (eg different filtering parameters in preprocessing)
    
    Input should be a directory with containing preprocessed h5ad files.

    Required annotations:
        - adata.obs["n_counts"], if model is cellassign (used as internal model parameter)

    Outputs a gzip compressed h5ad file for each file, with cell type annotations added to each cell.
    Output files are named {output_prefix}_{basename}_{i}.h5ad.

    Conditional annotations added to adata.obs: 
    - if model is z_score:
        [numpy.float32: "<cell_type>_score", str: "cell_type"] 
        (score for each cell type denoted in the marker file, chosen cell type for each cell)
    - if model is cellassign:
        [str: "cell_type"]
        (cell type prediction for each cell)

    Parameters:
        input_data_dir (str): path to directory containing h5ad files to annotate.
        use_ensembl_ids (bool): whether ensembl ids have been used for gene names.
        marker_file_path (str): path to file with marker genes for each cell type.
            This should be a json with cell type names as keys and lists of marker gene symbols as values.
        negative_marker_file_path (str, optional): path to file with negative marker genes for each cell type.
            This should be a json with cell type names as keys and lists of marker gene symbols as values.
            Only used if model is "z_score".
        model (Literal["z_score","cellassign"], optional): which model to use. Defaults to "cellassign".
        cutoff_unsure (float, optional): a value between 0 and 1, specifying how high the second highest 
            score can at most be relative to the highest to still annotate a well defined cell type
            Defaults to 0.8 (0.8 times highest score). Only used if model is "z_score".
        cutoff_other (float, optional): how many stdevs above / below the mean score the lowest cell type 
            score has to be to annotate a well defined cell type. Defaults to -0.2. Only used if model is "z_score".
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/cell_type_annotated. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "preprocessed".
        output_prefix (str, optional): prefix for output file names. Defaults to "cell_type_annotated".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have cell_type_annotated folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "cell_type_annotated"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated")
    output_temp_dir = os.path.join(TEMP_DIR, "cell_type_annotated")

    # assign output file list
    output_file_list = []

    # run script on each file in input_data_dir
    for file in os.listdir(input_data_dir):
        print("Annotating cell types for: " + file)
        temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Cell_type_annotation.py"), os.path.join(input_data_dir, file), output_temp_dir, [use_ensembl_ids, marker_file_path, negative_marker_file_path, model, cutoff_unsure, cutoff_other, verbose])
        
        # rename output file
        os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{file.removeprefix(input_prefix + "_")}"))
        temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{file.removeprefix(input_prefix + "_")}")

        # add output file to output_file_list
        output_file_list.append(temp_output_path)

        # if specified, permanently store a copy of the temporary output file
        if save_output == True:
            print(f"Saving {temp_output_path} to {output_storage_dir}")
            shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def aggregate_batches(
        input_data_dir: str,  
        max_obs_cancerous: int = None,
        max_obs_non_cancerous: int = None,
        save_output: bool = False, 
        input_prefix: str = "cell_type_annotated",
        output_prefix: str = "aggregated",
        verbose: bool = False) -> list[str]:
    """ 
    Runs batch_aggregation.py on input_data_dir to create one 
    aggregated h5ad file from all contained h5ad files of one cancer type (e.g. PDAC, reads from file name). Keeps all obs annotations.
    Tries to keep all present var annotations, but leaves out those that have non 
    unique var_name -> value mappings. 

    Input should be a directory with containing h5ad files that should be aggregated into one.

    Outputs a gzip compressed aggregated h5ad file.
    Output files are named {output_prefix}_{basename}.h5ad. (index is omitted because of aggregation)

    Annotations that are added to adata.obs: [str: "batch", str: "cancer_state"] (orginal batch of cell, whether cell is cancerous or not)

    Parameters:
        input_data_dir (str): path to directory containing h5ad files to annotate.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/cell_type_annotated. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "cell_type_annotated".
        output_prefix (str, optional): prefix for output file names. Defaults to "aggregated".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have cell_type_annotated folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "aggregated"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "aggregated"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "aggregated")
    output_temp_dir = os.path.join(TEMP_DIR, "aggregated")

    # assign output file list
    output_file_list = []

    # run script on each file in input_data_dir
    print("Aggregating files in " + os.path.basename(input_data_dir))
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "batch_aggregation.py"), input_data_dir, output_temp_dir, [max_obs_cancerous, max_obs_non_cancerous, input_prefix, output_prefix, verbose])
    
    # naming happens in subprocess as it relies on the isolated cancer type (extracted from filename)

    # add output files to output_file_list (here, the temp_output_path is a directory, so get its contents)
    for file in os.listdir(temp_output_path):
        output_file_list.append(os.path.join(temp_output_path, file))

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        for file in os.listdir(output_temp_dir):
            print(f"Saving {file} to {output_storage_dir}")
            shutil.copy(os.path.join(output_temp_dir, file), os.path.join(output_storage_dir, file))

    return output_file_list

def correct_batch_effects(
        input_data_file: str,
        max_considered_genes: int | Literal["all"] = 3000,
        save_output: bool = False,
        input_prefix: str = "aggregated",
        output_prefix: str = "batch_corrected",
        verbose: bool = False) -> list[str]:
    """
    Runs Batch_correction.py on an aggregated h5ad file to correct batch effects.

    One scvi and once scanvi model are trained once in the first run.
    The models are reused for all subsequent runs.
    To guarantee reproducibility, use the exact same input files for each run.

    Input should be an aggregated h5ad file.

    Outputs a gzip compressed corrected h5ad file with batch corrected representations.
    Output files are named {output_prefix}_{basename}.h5ad.
    
    Annotations added to adata.obsm: [pandas.DataFrame: "X_scVI_corrected", pandas.DataFrame: "X_scANVI_corrected"] (scVI and scANVI corrected normalized representations)

    Parameters:
        input_data_file (str): path to aggregated h5ad file to correct.
        max_considered_genes (int, optional): maximum number of highly variable genes to consider for model training. Defaults to 1000.
            if this is set to "all", scVI and scANVI will be run on the full geneset (might take a long time).
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/batch_corrected. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "aggregated".
        output_prefix (str, optional): prefix for output file names. Defaults to "batch_corrected".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "batch_corrected"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected")
    output_temp_dir = os.path.join(TEMP_DIR, "batch_corrected")
    
    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Correcting batch effects in {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Batch_correction.py"), input_data_file, output_temp_dir, [max_considered_genes, verbose])

    # rename output file
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def infer_CNVs(
        input_data_file: str, 
        reference_genome_path: str, 
        corrected_representation: str = None, 
        cell_type: str = None,
        save_output: bool = False,
        input_prefix: str = "batch_corrected",
        output_prefix: str = "CNV_inferred",
        verbose: bool = False) -> list[str]:
    """ 
    Infer copy number variations from any aggregated h5ad file.
    This works much better with with larger genesets, because it uses a sliding 
    window to pass over genes that are close on the chromosome.
    If cell_type is passed, it also reduce the adata to just that cell type.

    Input should be an h5ad file, batch corrected if corrected_representation is passed, cell type annotated if cell_type is passed.

    Requires the following annotations to be present:
        - if cell_type is passed, adata.obs["cell_type"] (from cell_type_annotation.py)
        - if corrected_representation is passed, adata.obsm[corrected_representation] (from batch_correction.py)
        - adata.obs["cancer_state"] (from batch_aggregation.py) (for root cell selection, and cancaer_state_inferred annotation)

    Outputs a gzip compressed corrected h5ad file with CNV annotations, reduced to cell type if cell_type is passed.
    Output files are named {output_prefix}_{basename}{suffix}.h5ad.
    where suffix is "_{cell_type}" if cell_type is passed

    Annotations added to adata.obs: [numpy.float64:"summed_cnvs", numpy.float64:"cnv_score", str:"cancer_state_inferred"] (number of CNVs in a cell, cnv score as computed by inferCNVpy.tl.cnv_score, whether a cell has has cnv_score > expected_cancer_percentage (expected_cancer_percentage is the percentage of cancer cells in adata.obs["cancer_state"],values are "non_cancerous" or "cancerous"))
    Annotations added to adata.obsm: [scipy.sparse._csr.csr_matrix[numpy.float64]:"{corrected_representation}_cnv", numpy.ndarray[numpy.float64]:"{corrected_representation}_gene_values_cnv"] 
    (smoothed and denoised gene expression along genomic locations (basically "cnvness of a gene"), number of copies per gene)

    Parameters:
        input_data_file (str): path to aggregated / batch corrected h5ad file.
        reference_genome_path (str): path to reference genome in gtf format. Eg. hg38.gtf.gz.
            Required to map var_names (which must be ensembl IDs, to allow unique mapping) to chromosomal coordinates.
        corrected_representation (str, optional): name of corrected representation on adata.obsm to use. 
            Defaults to None, meaning adata.X will be used.
        cell_type (str, optional): cell type (from adata.obs["cell_type"]) to infer CNVs for. Defaults to none, meaning cnvs are inferred for all cells.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/CNV. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "batch_corrected".
        output_prefix (str, optional): prefix for output file names. Defaults to "CNV_inferred".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "CNV"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "CNV"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "CNV")
    output_temp_dir = os.path.join(TEMP_DIR, "CNV")

    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Inferring CNVs from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "infer_CNV.py"), input_data_file, output_temp_dir, [reference_genome_path, corrected_representation, cell_type, verbose])

    # rename output file
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}{"_" + cell_type if cell_type else ""}.h5ad"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}{"_" + cell_type if cell_type else ""}.h5ad")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def run_scMF(
        input_data_file: str,
        layer: str = "X",
        cell_type = None,
        save_output: bool = False,
        input_prefix: str = "reduced",
        output_prefix: str = "scMF",
        verbose: bool = False) -> list[str]:
    
    """
    Run scMalignantFinder on the spcified layer in an h5ad file. 
    The script can also isolate ductal cells first, if specified.

    Input should be an h5ad file. Should be batch corrected.

    Requires the following annotations to be present:
        - if adata.var_names ) ensembl IDs, adata.var["gene_symbols"]
        - if layer is passed and != "X", adata.layers[layer] / adata.obsm[layer]
        - if cell_type is passed, adata.obs["cell_type"] (from cell_type_annotation.py)

    Outputs a gzip compressed h5ad file with malignancy state and malignancy probability annotations.
    Output files are named {output_prefix}_{basename}{suffix}.h5ad.
    where suffix is "_{cell_type}" if cell_type is passed

    Annotations added to adata.obs: [str:"cancer_state_inferred_scMF", numpy.float64: "malignancy_probability_scMF"] (whether a cell is cancerous or non_cancerous, probability a cell is cancerous)

    Parameters:
        input_data_file (str): path to aggregated / batch corrected h5ad file.
        layer (str, optional): layer / obsm matrix to run scMalignantFinder on. Defaults to "X", meaning adata.X will be used.
        cell_type (str, optional): cell type (from adata.obs["cell_type"]) to infer CNVs for. Defaults to none, meaning cnvs are inferred for all cells.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/CNV. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "batch_corrected".
        output_prefix (str, optional): prefix for output file names. Defaults to "scMF".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "scMF"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "scMF"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "scMF")
    output_temp_dir = os.path.join(TEMP_DIR, "scMF")

    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Running scMF on {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "sc_malignant_finder.py"), input_data_file, output_temp_dir, [layer, cell_type, verbose])

    # rename output file
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}{"_" + cell_type if cell_type and (cell_type not in os.path.basename(input_data_file)) else ""}.h5ad"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}{"_" + cell_type if cell_type and (cell_type not in os.path.basename(input_data_file)) else ""}.h5ad")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def get_phylogenetic_tree(
        input_data_file: str, 
        cnv_score_matrix: str,
        distance_metric: Literal["euclidean", "correlation"] = "correlation",
        n_clades: int = 30,
        grouping_metric: str = "cancer_state_inferred",
        transition_entropy_threshold: float = 0.8,
        save_output: bool = False,
        input_prefix: str = "scMF",
        output_prefix: str = "transition_clades",
        verbose: bool = False) -> list[str]:
    
    """
    Create a neighbor joining tree of correlation distances between CNV profiles of cells.
    Add non_cancerous / transitional / cancerous annotations to cells. 
    Cells are labelled normal / cancer if the entropy of the clade is below threshold based on dominant cell state in the calde.
    Some cells will not be assigned to a clade (due to the nature of the algorithm).
    These cells will be assigned to the virtual calde "-1" instead, and are annotated with "unassigned" 

    Input should be an h5ad file containing the cnv_score matrix and the obs annotation that will be used to color the leaves of the tree.


    Annotations that need to be present:
        - adata.obs[grouping_metric]
        - adata.obsm[cnv_score_matrix]

    Outputs newick tree file and h5ad file with new annoations. Returns directory these files are temprorarily saved in.
    Files saved permanently to OUTPUT_STORAGE_DIR/tree if save_output == True.
    
    Annotations added to adata.obs: [np.float: "cnv_clade", str:"cancer_state_inferred_tree"] (the clade a cell is in, on the phylogenetic tree (interger values from 0 to n_clades, or -1); a cell's inferred cancer state from the tree, includes "non_cancerous", "transitional", "cancerous", "unassigned")


    Parameters:
        input_data_file (str): path to h5ad file.
        cnv_score_matrix (str): name of cnv_score matrix in adata.obsm.
        distance_metric (Literal["euclidean", "correlation"], optional): metric used to compute distance between cells' cnv profiles for neighbor joining tree. Defaults to "correlation".
        n_clades (int, optional): number of clades to split tree into for transition state determination. Defaults to 30.
        grouping_metric (str, optional): adata.obs annotation used to determine cancer vs non cancer cells for transition state determination. Defaults to "cancer_state_inferred".
        transition_entropy_threshold (float, optional): Clades with entropy (mixing of cancerous vs non_cancerous cells according to grouping_metric) above this value are considered transitional. Defaults to 0.8.
        save_output (bool, optional): whether to save the plots permanently. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "reduced".
        output_prefix (str, optional): prefix for output file names. Defaults to "transition_clades".
        verbose (bool, optional): whether to print verbose output. Defaults to False.

    Returns:
        list[str]: list of paths to output files. One h5ad file and one newick tree file.
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have relevant folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "tree"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "tree"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "tree")
    output_temp_dir = os.path.join(TEMP_DIR, "tree")

    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Generating phylogenetic tree: {input_data_file}")
    hf.execute_subprocess(os.path.join(SCRIPT_DIR, "phylogenetic_tree.py"), input_data_file, output_temp_dir, [cnv_score_matrix, distance_metric, n_clades, grouping_metric, transition_entropy_threshold, verbose])

    # rename output file (this is the h5ad file) (tree naming happens in subprocess) 
    os.rename(os.path.join(output_temp_dir, os.path.basename(input_data_file)), os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}.h5ad"))

    # add output files to output_file_list
    for output_file in os.listdir(output_temp_dir):
        output_file_list.append(os.path.join(output_temp_dir, output_file))

    # if specified, permanently store a copy of the temporary output files from the output_temp_dir
    if save_output == True:
        for file in os.listdir(output_temp_dir):
            shutil.copy(os.path.join(output_temp_dir, file), os.path.join(output_storage_dir, file))

    return output_file_list # contains h5ad file and newick tree file paths


def isolate_and_HVGs(
        input_data_file: str, 
        main_layer: str = None,
        add_log1p: bool = False,
        max_considered_genes: int | Literal["all"] = 3000,
        genes_to_keep: list[str] = [],
        isolation_dict: dict[str, list] = {},
        preservation_dict: dict[str, list] = {},
        save_output: bool = False,
        input_prefix: str = "transition_clades",
        output_prefix: str = "isolated",
        verbose: bool = False) -> list[str]:
    """ 
    Move main layer to adata.X, purge everything else. Optionally add log1p layer of main layer.
    Optionally isolate cells according to isolation_dict.
    Optionally compute HVGs (are computed on log1p layer, if log1p layer is added). 
    In this order.

    Input should be an h5ad file.

    Requires the following annotations to be present: 
        - the main_layer (if main_layer is not None) in adata.X, adata.layers or adata.obsm
        - all conditions in isolation_dict must have a column in adata.obs with that has at least one entry of the type to isolate (e.g. cancer_state: ["normal", "transitional"] -> requires adata.obs["cancer_state"].unqiue() == ["normal", "transitional", ...])
        - all conditions in preservation_dict must have a column in adata.obs with that has at least one entry of the type to isolate (e.g. cancer_state: ["normal", "transitional"] -> requires adata.obs["cancer_state"].unqiue() == ["normal", "transitional", ...])
        - all values in genes_to_keep must exist in adata.var_names, adata.var["gene_symbols"] or adata.var["gene_ids"]

    Outputs a gzip compressed h5ad file with, optionally, reduced layers / obsms and HVGs.
    Output files are named {output_prefix}_{basename}{suffix}.h5ad.
    where suffix may contain "_HVG" if HVGs are selected, "_X_is_{main_layer}" if main layer is given and "_{isolation_dict.key}_is_{isolation_dict[key]}" for each key in isolation_dict

    Annotations added to adata.layers:
        - if add_log1p == True: adata.layers["log1p"] = log1p(main_layer)

    Parameters:
        input_data_file (str): path to h5ad file.
        main_layer (str, optional): layer / obsm that will be moved to adata.X. Defaults to None, meaning adata.X will remain as is.
        add_log1p (bool, optional): whether to add a new layer with log1p(main_layer). Defaults to False.
        max_considered_genes (int, optional): maximum number of HVGs to consider. Defaults to 3000. If this is set to "all", no HVGs will be selected.
        genes_to_keep (list[str], optional): these genes will be kept regardless of their HVG status. Defaults to [].
        isolation_dict (dict[str, list], optional): dictionary of conditions to isolate, each key is a condition and each value is a list of the types to isolate. Defaults to {}
        preservation_dict (dict[str, list], optional): dictionary of conditions where cells that have them will be kept regardless of isolation_dict (eg. {"cnv_clade" : 28} -> all cells from cnv clade 28 will be kept). Defaults to {}.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/reduced. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "batch_corrected".
        output_prefix (str, optional): prefix for output file names. Defaults to "CNV_inferred".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Returns:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "isolated"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "isolated"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "isolated")
    output_temp_dir = os.path.join(TEMP_DIR, "isolated")

    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Reducing {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "matrix_isolation_HVGs.py"), input_data_file, output_temp_dir, [main_layer, add_log1p, max_considered_genes, genes_to_keep, isolation_dict, preservation_dict, verbose])

    # rename output file
    HVG_suffix = "_HVG" if max_considered_genes != "all" else ""
    X_is_suffix = f"_X_is_{main_layer}" if main_layer is not None else ""
    isolation_suffix = ("_" + "_".join(f"{key}_is_{value}" for key, value in isolation_dict.items())) if isolation_dict != {} else ""
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}{HVG_suffix}{X_is_suffix}{isolation_suffix}.h5ad"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_").removesuffix(".h5ad")}{HVG_suffix}{X_is_suffix}{isolation_suffix}.h5ad")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list

def infer_pseudotime(
        input_data_file: str, 
        origin_clade: int,
        flavor: Literal["monocle", "dpt"] = "monocle",
        layer: str = "log1p",
        smoothe_expression: bool = True,
        find_switches: bool = True,
        bic_threshold: float = 10,
        mean_threshold: float = 0.5,
        save_output: bool = False,
        input_prefix: str = "isolated",
        output_prefix: str = "pseudotime_inferred",
        verbose: bool = False) -> list[str]:
    """ 
    Infer pseudotime from an h5ad file. Optionally smoothe expression per gene along pseuodtime with sliding window smoothing.

    Input should be a an h5ad file. Ideally with transition clades isolated.

    Requires the following annotations to be present:
        - if origin_clade is passed, adata.obs["cnv_clades"] (from phylogenetic_tree.py)
        - adata.layers[layer] / adata.obsm[layer] / adata.X; expression matrix that should be used for pseudotime
          (log1p is recommended)

    Outputs a gzip compressed h5ad with pseudotime annotated for each cell. 
    Output files are named {output_prefix}_{basename}.h5ad.

    Annotations added to adata.obs:
        - if flavor == "dpt", adata.obs['dpt_pseudotime']: pandas.Series (dtype float) (pseudotime value for each cell)
        - if flavor == "monocle", adata.obs['monocle_pseudotime']: pandas.Series (dtype float) (pseudotime value for each cell)
        
    Annotations added to adata.var:
        - if find_switches == True, adata.var['is_switch']: bool (wether the gene is a switch or not)
        - if find_switches == True, adata.var['bic_k1']: float (bic for 1 component GMM)
        - if find_switches == True, adata.var['bic_k2']: float (bic for 2 component GMM)

    Parameters:
        input_data_file (str): path to aggregated / batch corrected h5ad file.
        origin_clade (int, optional): root clade for pseudotime inference (this should be the most normal clade, lowest cnv / whatver fitting metric you use).
        flavor (Literal["monocle", "dpt"], optional): pseudotime inference method. Defaults to "monocle".
        layer (str, optional): layer / obsm / adata.X to use for pseudotime inference. Defaults to "log1p".
        smoothe_expression (bool, optional): whether to smoothe expression matrix along pseudotime. Defaults to True.
        find_switches (bool, optional): whether to isolate genes that switch expression along pseudotime (ON / OFF). Defaults to True.
        bic_threshold (float, optional): genes are considered switching genes if the bic for a 2 component GMM is at most the bic for the 1 component GMM - bic_threshold. Set higher for a stricter selection. Defaults to 10.
        mean_threshold (float, optional): switching genes identified by the GMM method are only get is_switch = True if the mean expression of the higher component is at least mean(lower component) + mean(lower component)*mean_threshold. Defaults to 0.5.
        save_output (bool, optional): whether to save output files permanently to OUTPUT_STORAGE_DIR/CNV. Defaults to False.
        input_prefix (str, optional): prefix of input file names, must match or will cause error. Defaults to "batch_corrected".
        output_prefix (str, optional): prefix for output file names. Defaults to "CNV_inferred".
        verbose (bool, optional): whether to print verbose output from subprocess. Defaults to False.

    Return:
        list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "pseudotime"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "pseudotime"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "pseudotime")
    output_temp_dir = os.path.join(TEMP_DIR, "pseudotime")

    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Inferring pseudotime from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "pseudotime_inference.py"), input_data_file, output_temp_dir, [origin_clade, flavor, layer, smoothe_expression, find_switches, bic_threshold, mean_threshold, verbose])

    # rename output file
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def infer_GRN_edges(
        input_data_file: str,
        layer = "log1p",
        n_nodes = None,
        tf_list_file: str = None, 
        convergence_threshold: float = 0.05,
        top_n_regulators: int = 5,
        min_runs: int = 20,
        min_stability: int = 0.5,
        verbose: bool = False,
        input_prefix: str = "pseudotime_inferred",
        output_prefix: str = "GRN_edges",
        save_output: bool = False) -> list[str]:
    """ 
    Infer the edges of a GRN from a given h5ad file by bootstrapping grnboost2.
    A transcription factor list can be passed to limit the space of potential regulators.

    Input should be an h5ad file, ideally with expression smoothed along pseudotime and switch genes isolated.
    
    Requires the following annotations to be present:
        - adata.layers[layer]; expression matrix that should be used to infer edges from.
          (log1p wtih smoothed expression is recommended)

    Outputs a json file where keys = target gene, values = list of most stable regulators

    This script adds no annotations to adata, as it does not return an h5ad file.

    Parameters
    ----------
    input_data_file : str
        path to h5ad file.
    layer : str, optional
        layer in adata.layers that should be used for edge inference. Defaults to "log1p".
        Will use adata.X if layer is "X".
    n_nodes : int, optional
        maximum number of nodes to consider for GRN. Defaults to None meaning all genes are considered as nodes.
        (this is primarily used for debugging as it randomly chooses genes from adata.var_names)
    tf_list_file : str, optional
        path to transcription factor list file (txt). Defaults to None, meaning all genes will be considered as potential regulators.
    convergence_threshold : float, optional
        convergence limit; stops when average edge stability (stdev of occurance across runs) relatively changes by less than this value from one run to next. Defaults to 0.05 = 5%.
    top_n_regulators : int, optional
        maximum number of regulators to consider for each target gene. Defaults to 5.
    min_runs : int, optional
        minimum amount of runs before converging (it could otherwise happen that stdevs randomly fall close to eachother and premture false convergence occurs). Defaults to 20.
    min_stability : float, optional
        putative regulators need to occur at least in this percentage of rens eg 0.5 = 50% of runs. Defaults to 0.5.
    verbose: bool, optional
        whether to print verbose output from subprocess.
    input_prefix : str, optional
        prefix of input file names, must match or will cause error. Defaults to "pseudotime".
    output_prefix : str, optional
        prefix for output file names. Defaults to "GRN_edges".
    save_output : bool, optional
        whether to save output files permanently to OUTPUT_STORAGE_DIR/GRN_edges. Defaults to False.
        
    Returns
    -------
    list[str]: list of paths to output files
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "GRN_edges"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "GRN_edges"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "GRN_edges")
    output_temp_dir = os.path.join(TEMP_DIR, "GRN_edges")

    # assign output file list
    output_file_list = []

    # run script and assign path to temporary output file
    print(f"Inferring GRN edges from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "GRN_edge_inference.py"), input_data_file, output_temp_dir, [layer, n_nodes, tf_list_file, convergence_threshold, top_n_regulators, min_runs, min_stability, input_prefix, output_prefix, verbose])

    # naming happens in subprocess
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list


def infer_GRN_rules(input_data_file: str, edge_set_file: str, verbose: bool = False):
    """ 
    Infer the rules of a GRN from a given aggregated / batch corrected h5ad file
    and a JSON with inferred edges for each target gene. Outputs a JSON with inferred rule
    for each target gene.
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "GRN_rules"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "GRN_rules"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "GRN_rules")
    output_temp_dir = os.path.join(TEMP_DIR, "GRN_rules")

    # run script and assign path to temporary output file
    print(f"Inferring GRN from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "GRN_rule_inference.py"), input_data_file, output_temp_dir, [verbose, edge_set_file])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["GRN_rule_inference.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


def simulate_GRN(input_data_file: str, verbose: bool = False):
    """ 
    Simulate an existing GRN and do a modular analysis thereof. Input is a txt with
    a lost of boolean rules in BNET format. 
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "GRN_analysis"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "GRN_analysis"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "GRN_analysis")
    output_temp_dir = os.path.join(TEMP_DIR, "GRN_analysis")

    # run script and assign path to temporary output file
    print(f"Inferring GRN from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "GRN_simulation.py"), input_data_file, output_temp_dir, [verbose])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["GRN_simulation.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))            


# main loop
if __name__ == "__main__": # ensures this code runs only when this script is executed directly, not when imported
    
    # --- register signal handlers ---
    import signal
    signal.signal(signal.SIGINT, request_stop)   # Ctrl+C
    signal.signal(signal.SIGTERM, request_stop)  # kill (default TERM)
    try:
        signal.signal(signal.SIGHUP, request_stop)   # terminal close (Linux/macOS)
    except AttributeError:
        pass  # SIGHUP not available on Windows

    # --- main loop ---

    try:
        # full pipeline run (adjust parameters as needed)
        r"""use_ensembl_ids = True
        mode = choose_pipeline_mode(RAW_DATA_DIRS[0])
        for data_dir in RAW_DATA_DIRS:
            preprocess_data(data_dir, mode, use_ensembl_ids=use_ensembl_ids, save_output=True, verbose=True, filtering_params=FilteringParameters(min_n_cells_percentage=0, max_n_MADs=None, max_mito_percentage=0.10))
        annotate_cell_types(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed"), use_ensembl_ids, os.path.join(AUX_DATA_DIR, "annotations", "marker_genes.json"), model="cellassign", verbose=True, save_output=True)
        output_path_list = aggregate_batches(os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated"), save_output=True, verbose=True)
        output_path_list = correct_batch_effects(output_path_list[0], max_considered_genes="all", save_output=True, verbose=True)
        output_path_list = infer_CNVs(output_path_list[0], corrected_representation="X_scANVI_corrected", reference_genome_path=os.path.join(AUX_DATA_DIR, "annotations", "gencode.v49.annotation.gtf.gz"), cell_type="ductal_cell", save_output=True, verbose=True)
        output_path_list = reduce_data(output_path_list[0], input_prefix="CNV_inferred", layers_to_remove=["X", "X_scANVI_corrected_gene_values_cnv", "X_scVI_corrected"], save_output=True, verbose=True)
        output_path_list = run_scMF(output_path_list[0], layer="X_scANVI_corrected", save_output=True, verbose=True)
        output_path_list = get_phylogenetic_tree(output_path_list[0], cnv_score_matrix="X_scANVI_corrected_cnv", distance_metric="euclidean", n_clades=30, grouping_metric="cancer_state_inferred_scMF", transition_entropy_threshold=0.8, save_output=True, verbose=True)
        tree_file = None
        h5ad_file = None
        for file in output_path_list:
            if ".nwk" in file:
                tree_file = file
            elif ".h5ad" in file:
                h5ad_file = file
        #plots for cnv cluster validation
        for projection in ["UMAP", "PCA"]:
            cluster_and_plot(["projections"], input_data_file=h5ad_file, obs_annotations=["cancer_state", "cancer_state_inferred", "cancer_state_inferred_scMF", "cancer_state_inferred_tree", "cnv_score", "cnv_clade"], layers=["X_scANVI_corrected_cnv"], projection=projection, output_storage_subdir="CNV_matrix_check", save_output=True, verbose=True, show=False)                
        # tree plot (CANNOT RUN ON CLUSTER, BCS PLOTLY NEEDS MANUAL SAVING)
        # cluster_and_plot(["phylogenetic_tree"], input_data_file=h5ad_file, tree_file=tree_file, output_storage_subdir="tree", save_output=True, verbose=True, show=False)
        output_path_list = isolate_and_HVGs(h5ad_file, main_layer="X_scANVI_corrected", add_log1p=True, max_considered_genes=3000, isolation_dict={"cancer_state_inferred_tree":["transitional"]}, save_output=True, verbose=True)
        # plots for clade selection
        for projection in ["UMAP", "PCA"]:
            cluster_and_plot(["projections"], input_data_file=output_path_list[0], obs_annotations=["cancer_state", "cancer_state_inferred", "cancer_state_inferred_tree", "cnv_score", "cnv_clade"], layers=["log1p"], projection=projection, output_storage_subdir="clade_selection", save_output=True, verbose=True, show=False)
        """

        # subsample_cells(file_path=os.path.join(OUTPUT_STORAGE_DIR, "scMF", "scMF_run0_PDAC_ductal_cell.h5ad"), fraction=0.6, n_samples=3, output_dir=os.path.join(OUTPUT_STORAGE_DIR, "subsampled"))
        
        # identify clades to keep for pseudotime
        r"""        
        for run in os.listdir(os.path.join(OUTPUT_STORAGE_DIR, "subsampled")):
            output_path_list = get_phylogenetic_tree(os.path.join(OUTPUT_STORAGE_DIR, "subsampled", run), cnv_score_matrix="X_scANVI_corrected_cnv", distance_metric="euclidean", n_clades=30, grouping_metric="cancer_state_inferred", transition_entropy_threshold=0.8, save_output=True, verbose=True)
            for file in output_path_list:
                if ".nwk" in file:
                    tree_file = file
                elif ".h5ad" in file:
                    h5ad_file = file
            output_path_list = isolate_and_HVGs(h5ad_file, main_layer="X_scANVI_corrected", add_log1p=True, max_considered_genes=3000, isolation_dict={}, save_output=True, verbose=True)
            cluster_and_plot(["projections"], input_data_file=output_path_list[0], obs_annotations=["cancer_state", "cancer_state_inferred", "cancer_state_inferred_tree", "cancer_state_inferred_scMF" "cnv_score", "cnv_clade"], layers=["log1p"], projection="UMAP", output_storage_subdir="clade_selection", save_output=True, verbose=True, show=False)
        """

        # tree output -> grn edges
        r"""
        for run, clades in zip(os.listdir(os.path.join(OUTPUT_STORAGE_DIR, "tree")), [[18,9], [21,26], [26,19]]):
            if os.path.isdir(os.path.join(OUTPUT_STORAGE_DIR, "tree", run)):
                continue
            output_path_list = isolate_and_HVGs(os.path.join(OUTPUT_STORAGE_DIR, "tree", run), main_layer="X_scANVI_corrected", add_log1p=True, max_considered_genes=3000, isolation_dict={"cancer_state_inferred_tree":["transitional"]}, preservation_dict={"cnv_clade":clades}, save_output=True, verbose=True)
            output_path_list = infer_pseudotime(input_data_file=output_path_list[0], origin_clade=clades[0], save_output=True, verbose=True)
            infer_GRN_edges(input_data_file=output_path_list[0], tf_list_file=os.path.join(AUX_DATA_DIR, "annotations", "tf_symbols_list.txt"), verbose=True, save_output=True)
        """

        # verify ts consistency under different conditions
        r"""
        for file in os.listdir(os.path.join(OUTPUT_STORAGE_DIR, "subsampled")):
            output_path_list = get_phylogenetic_tree(os.path.join(OUTPUT_STORAGE_DIR, "subsampled", file), cnv_score_matrix="X_scANVI_corrected_cnv", distance_metric="euclidean", n_clades=30, grouping_metric="cancer_state_inferred_scMF", transition_entropy_threshold=0.8, save_output=True, verbose=True, output_prefix="transition_clades_scMF")
            for file in output_path_list:
                if ".nwk" in file:
                    tree_file = file
                elif ".h5ad" in file:
                    h5ad_file = file
            output_path_list = isolate_and_HVGs(h5ad_file, main_layer="X_scANVI_corrected", add_log1p=True, max_considered_genes=3000, isolation_dict={"cancer_state_inferred_tree":["transitional"]}, preservation_dict={}, save_output=True, verbose=True)
        """

        # check BC condition's effects on TS consistency
        r"""
        #output_path_list = correct_batch_effects(os.path.join(OUTPUT_STORAGE_DIR, "aggregated", "aggregated_run0_PDAC.h5ad"), max_considered_genes="all", save_output=False, verbose=True, output_prefix="batch_corrected_seed_test")
        #output_path_list = infer_CNVs(output_path_list[0], corrected_representation="X_scANVI_corrected", reference_genome_path=os.path.join(AUX_DATA_DIR, "annotations", "gencode.v49.annotation.gtf.gz"), cell_type="ductal_cell", save_output=False, verbose=True, input_prefix="batch_corrected_seed_test", output_prefix="CNV_inferred_seed_test")
        #output_path_list = reduce_data(output_path_list[0], layers_to_remove=["X", "X_scANVI_corrected_gene_values_cnv", "X_scVI_corrected"], save_output=False, verbose=True, input_prefix="CNV_inferred_seed_test", output_prefix="reduced_seed_test")
        output_path_list = get_phylogenetic_tree(os.path.join(OUTPUT_STORAGE_DIR, "CNV", "CNV_inferred_run0_PDAC_ductal_cell.h5ad"), cnv_score_matrix="X_scANVI_corrected_cnv", distance_metric="euclidean", n_clades=30, grouping_metric="cancer_state_inferred", transition_entropy_threshold=0.8, save_output=True, verbose=True, input_prefix="reduced_seed_test", output_prefix="transition_clades_no_changes")
        for file in output_path_list:
            if ".nwk" in file:
                tree_file = file
            elif ".h5ad" in file:
                h5ad_file = file
        output_path_list = isolate_and_HVGs(h5ad_file, main_layer="X_scANVI_corrected", add_log1p=True, max_considered_genes=3000, isolation_dict={"cancer_state_inferred_tree":["transitional"]}, preservation_dict={}, save_output=True, verbose=True)
        """

        for file in os.listdir(os.path.join(OUTPUT_STORAGE_DIR, "subsampled", "jaccard_convergency_check")):
            output_path_list = get_phylogenetic_tree(os.path.join(OUTPUT_STORAGE_DIR, "subsampled", "jaccard_convergency_check", file), cnv_score_matrix="X_scANVI_corrected_cnv", distance_metric="euclidean", n_clades=30, grouping_metric="cancer_state_inferred", transition_entropy_threshold=0.8, save_output=False, verbose=False, output_prefix="transition_clades_jaccard_convergence")
            for file in output_path_list:
                if ".nwk" in file:
                    tree_file = file 
                elif ".h5ad" in file:
                    h5ad_file = file
            output_path_list = isolate_and_HVGs(h5ad_file, main_layer="X_scANVI_corrected", add_log1p=True, max_considered_genes=3000, isolation_dict={"cancer_state_inferred_tree":["transitional"]}, preservation_dict={}, save_output=True, verbose=True)


        purge_tempfiles()
        sys.exit(0)
    except Exception:
        purge_tempfiles()
        raise # re-raise the exception to see the traceback and error message


