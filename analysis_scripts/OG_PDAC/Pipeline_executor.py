# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients
# change RAW_DATA_DIR and outcome_storage for each run to reflect new input and what outcomes should be permanenlty saved

import subprocess # needed for running other scripts
import os # needed for file and directory operations
from enum import Enum
import shutil # needed for file storage operations
import tempfile # needed for temporary file operations
import sys # needed to exit the program
from dataclasses import dataclass
import helper_functions as hf



#check if the required directories exist, if not create them
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data")):  # check if data directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data"))  # create data directory if it does not exist
    print("Created Data directory. Please add subdirectories for each sample containing the mtx and tsv outputs of the 10x genomics pipeline")  # inform user to add subdirectories with required files
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "output_storage")):  # check if storage directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "output_storage"))  # create storage directory if it does not exist
    print("Created output_storage directory. Intermediate output files will be saved here, if specified in OUTCOME_STORAGE.")  # inform user that intermediate outputs can be saved here
if not os.path.exists(os.path.join(tempfile.gettempdir(),"python")):
    os.makedirs(os.path.join(tempfile.gettempdir(),"python"))
    print("created \"python\" directory in appdata/local/temp to store temporary pipeline files. These files will be deleted when purgetempfiles() is called.")



# path constants
SCRIPT_DIR = os.path.dirname(__file__)  # directory where this script is located
# list of directories (see choose_pipeline_mode for valid structures for each entry)
RAW_DATA_DIRS = [    
                os.path.join(SCRIPT_DIR, "..", "..", "Data","pretraining", "cancerSCEM", "breast_cancer_cancerous")
                ]
OUTPUT_STORAGE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data", "output_storage")  # directory for optional permanent storage of indermediate subprocess outputs
TEMP_DIR = os.path.join(tempfile.gettempdir(),"python") # directory for storage of temporary pipeline files


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


@dataclass
class FilteringParameters:
    """
    Configuration for filtering low-quality cells and genes from single-cell data.

    Attributes:
        min_n_genes_percentile: 
            Cells with fewer expressed genes than this percentile are removed.
        min_n_cells_percentage: 
            Genes expressed in fewer than this fraction of cells are removed.
        min_n_UMIs_percentile:
            Cells with fewer UMI counts than this percentile are removed.
        max_n_MADs:
            Cells with mitochondrial percentages greater than (median + this many MADs) are removed.
        expected_doublet_percentage:
            Fraction of cells expected to be doublets, used by Scrublet. 
            See the sequencing device manufacturer’s recommendations.
    """

    # default values chosen based on manual inspection of plots of OG_PDAC data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966#:~:text=Summary%20Pancreatic%20ductal%20adenocarcinoma%20,plot%20to%20predict%20the%20overall
    min_n_genes_percentile: int = 10
    min_n_cells_percentage: float = 0.01
    min_n_UMIs_percentile: int = 10
    max_n_MADs: int = 1
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
            shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

        i += 1

    return output_file_list


def annotate_cell_types(
        input_data_dir: str, 
        marker_file_path: str, 
        cutoff_unsure: float = 0.8, 
        cutoff_other: float = -0.2, 
        use_ensembl_ids: bool = True, 
        save_output: bool = False, 
        input_prefix: str = "preprocessed",
        output_prefix: str = "cell_type_annotated",
        verbose: bool = False) -> list[str]:
    """ 
    Loops through input_data_dir and runs Cell_type_annotation.py on each contained h5ad file.

    Input should be a directory with containing preprocessed h5ad files.

    Outputs a gzip compressed h5ad file for each file, with cell type annotations added to each cell.
    Output files are named {output_prefix}_{basename}.h5ad.

    Annotations added to adata.obs: [numpy.float32: "<cell_type>_score", str: "cell_type"] 
    (score for each cell type denoted in the marker file, chosen cell type for each cell)

    Parameters:
        input_data_dir (str): path to directory containing h5ad files to annotate.
        marker_file_path (str): path to file with marker genes for each cell type.
            This should be a json with cell type names as keys and lists of marker gene symbols as values.
        cutoff_unsure (float, optional): a value between 0 and 1, specifying how high the second highest 
            score can at most be relative to the highest to still annotate a well defined cell type
            Defaults to 0.8 (0.8 times highest score).
        cutoff_other (float, optional): how many stdevs above / below the mean score the lowest cell type 
            score has to be to annotate a well defined cell type. Defaults to -0.2.
        use_ensembl_ids (bool, optional): whether ensembl ids have been used for gene names. Defaults to True.
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
        temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Cell_type_annotation.py"), os.path.join(input_data_dir, file), output_temp_dir, [marker_file_path, cutoff_unsure, cutoff_other, use_ensembl_ids, verbose])
        
        # rename output file
        os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{file.removeprefix(input_prefix + "_")}"))
        temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{file.removeprefix(input_prefix + "_")}")

        # add output file to output_file_list
        output_file_list.append(temp_output_path)

        # if specified, permanently store a copy of the temporary output file
        if save_output == True:
            shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

    return output_file_list

def correct_batch_effects(input_data_dir: str):
    """ WIP; DOES NOT WORK
    Run Batch_correction.py on a given directory of preprocessed h5ad files.
    Creates necessary directories if not present. Saves output permanenlty if specified."""

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "batch_corrected"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected")
    output_temp_dir = os.path.join(TEMP_DIR, "batch_corrected")

    # run script and assign path to temporary output file
    print(f"Correcting batch effects in {input_data_dir}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Batch_correction.py"), input_data_dir, output_temp_dir)

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["Batch_correction.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


def infer_CNVs(input_data_file: str, reference_genome_path: str, corrected_representation: str = None, verbose: bool = False):
    """ 
    Infer copy number variations from a given aggregated / batch corrected h5ad file.
    Requires to specify a reference genome, e.g. hg38 to map ensembl IDs (should be adata.var_names)
    to chromosomal coordinates. The reference genome should be in gtf format (gz compression is supported).
    Optionally pass a corrected representation in adata.ombsm (e.g. "X_scVI" or "X_scANVI").
    If not specified, the raw representation is used.
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "CNV"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "CNV"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "CNV")
    output_temp_dir = os.path.join(TEMP_DIR, "CNV")

    # run script and assign path to temporary output file
    print(f"Inferring GRN from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "infer_CNV.py"), input_data_file, output_temp_dir, [reference_genome_path, corrected_representation, verbose])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["infer_CNV.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


def prepare_for_pseudotime(input_data_dir: str):
    script_path = os.path.join(SCRIPT_DIR, "prepare_for_pseudotime.py")
    script_name = os.path.basename(script_path).removesuffix(".py")

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have relevant folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, script_name), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, script_name), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, script_name)
    output_temp_dir = os.path.join(TEMP_DIR, script_name)

    # run script and assign path to temporary output file
    print(f"Prepping for pseudotime: {input_data_dir}")
    hf.execute_subprocess(script_path, input_data_dir, output_temp_dir)

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE[script_name + ".py"] == True:
        for file in os.listdir(output_temp_dir):
            shutil.copy(os.path.join(output_temp_dir, file), os.path.join(output_storage_dir, os.path.basename(file)))

    return None


def infer_pseudotime(input_data_file: str, verbose: bool = False):
    """ 
    Infer pseudotime from a given aggregated / batch corrected h5ad file.
    Adds adata.obs["dpt_pseudotime"]
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "pseudotime"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "pseudotime"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "pseudotime")
    output_temp_dir = os.path.join(TEMP_DIR, "pseudotime")

    # run script and assign path to temporary output file
    print(f"Inferring GRN from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "pseudotime_inference.py"), input_data_file, output_temp_dir, [verbose])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["pseudotime_inference.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


def cluster_and_plot(input_data_dir: str, annotations: list = None, embedding: str = None, verbose: bool = False):
    if annotations is None:
        raise ValueError("annotations must be a list of strings describing adata.obs[...] column names")

    script_path = os.path.join(SCRIPT_DIR, "Clustering.py")
    script_name = os.path.basename(script_path).removesuffix(".py")

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have relevant folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, script_name), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, script_name), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, script_name)
    output_temp_dir = os.path.join(TEMP_DIR, script_name)

    # run script and assign path to temporary output file
    print(f"Clustering and plotting: {input_data_dir}")
    hf.execute_subprocess(script_path, input_data_dir, output_temp_dir, [annotations, embedding, verbose])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE[script_name + ".py"] == True:
        for file in os.listdir(output_temp_dir):
            shutil.copy(os.path.join(output_temp_dir, file), os.path.join(output_storage_dir, file))

    return None


def infer_GRN_edges(input_data_file: str, n_nodes = None, verbose: bool = False):
    """ 
    Infer the edges of a GRN from a given aggregated / batch corrected h5ad file.
    Outputs a JSON with inferred edges for each target gene as list.
    Set n_nodes to limit the number of nodes in the output GRN, by default, all
    genes are used.
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "GRN_edges"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "GRN_edges"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "GRN_edges")
    output_temp_dir = os.path.join(TEMP_DIR, "GRN_edges")

    # run script and assign path to temporary output file
    print(f"Inferring GRN from {input_data_file}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "GRN_edge_inference.py"), input_data_file, output_temp_dir, [verbose, n_nodes])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["GRN_edge_inference.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


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


def compute_variance():
    """Loops through the epithelial isolated h5ad files in temp/epithelial_isolated and runs Variance.py on each."""

    ADJ_paths = []  # list to store paths to tumor adjacent tissue files
    PDAC_paths = []  # list to store paths to tumor tissue files

    for file in os.listdir(os.path.join(OUTPUT_STORAGE_DIR, "epithelial_isolated")):  # iterate over files in the temp directory where epithelial isolated files are stored
        if file.endswith(".h5ad") and "ADJ" in file:  # check if the file is a h5ad file and tumor adjacent tissue
            ADJ_paths.append(os.path.join(OUTPUT_STORAGE_DIR, "epithelial_isolated", file))  # add the file path to the list
        elif file.endswith(".h5ad") and "PDAC" in file:  # check if the file is a h5ad file and tumor tissue
            PDAC_paths.append(os.path.join(OUTPUT_STORAGE_DIR, "epithelial_isolated", file))
        else:
            print(f"Variance: Skipping file {file} as it is not a .h5ad.gz file.")

    ADJ_paths = ",".join(ADJ_paths)  # join the list of paths into a single string separated by commas (so it can be passed as a command line argument)
    PDAC_paths = ",".join(PDAC_paths)

    subprocess.run(["python", os.path.join(SCRIPT_DIR, "Variance.py"), ADJ_paths, PDAC_paths], check=True)
            


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
    use_ensebml_ids = True # define whether to use ensembl ids, used for entire pipeline to avoid mismatches

    try:
        """for raw_data_dir in RAW_DATA_DIRS:
            mode = choose_pipeline_mode(raw_data_dir)
            preprocess_data(raw_data_dir, mode, use_ensembl_ids=use_ensebml_ids, save_output=False, verbose=True)
            """
        
        annotate_cell_types(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed"), marker_file_path=r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\marker_genes.json", use_ensembl_ids=use_ensebml_ids, save_output=True, verbose=True)
        # correct_batch_effects(os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated"))
        # infer_CNVs(os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected", "batch_corrected_HVG_PDAC.h5ad"), r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\gencode.v49.annotation.gtf.gz", corrected_representation="X_scANVI_corrected", verbose=True)
        # for file in ["annotated_PDAC_cancerous_4.h5ad", "annotated_PDAC_non_cancerous_2.h5ad"]:
        #     cluster_and_plot(os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated", file), ["cell_type"], verbose=True)
        
        purge_tempfiles()
        sys.exit(0)
    except Exception:
        purge_tempfiles()
        raise # re-raise the exception to see the traceback and error message