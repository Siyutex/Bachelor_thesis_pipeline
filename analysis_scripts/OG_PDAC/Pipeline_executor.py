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
from typing import Literal



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
                os.path.join(SCRIPT_DIR, "..", "..", "Data","OG_data", "NCBI","PDAC_cancerous"),
                os.path.join(SCRIPT_DIR, "..", "..", "Data","OG_data","NCBI","PDAC_non_cancerous")
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
            print(f"Saving {temp_output_path} to {output_storage_dir}")
            shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

        i += 1

    return output_file_list


def annotate_cell_types(
        input_data_dir: str, 
        use_ensembl_ids: bool, 
        marker_file_path: str, 
        negative_marker_file_path: str = None,
        model: Literal["z_score","cellassign"] = "z_score",
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
        (ratio of a cell's UMI counts to average cell's counts, cell type prediction for each cell)

    Parameters:
        input_data_dir (str): path to directory containing h5ad files to annotate.
        use_ensembl_ids (bool): whether ensembl ids have been used for gene names.
        marker_file_path (str): path to file with marker genes for each cell type.
            This should be a json with cell type names as keys and lists of marker gene symbols as values.
        negative_marker_file_path (str, optional): path to file with negative marker genes for each cell type.
            This should be a json with cell type names as keys and lists of marker gene symbols as values.
        model (Literal["z_score","cellassign"], optional): which model to use. Defaults to "z_score".
        cutoff_unsure (float, optional): a value between 0 and 1, specifying how high the second highest 
            score can at most be relative to the highest to still annotate a well defined cell type
            Defaults to 0.8 (0.8 times highest score).
        cutoff_other (float, optional): how many stdevs above / below the mean score the lowest cell type 
            score has to be to annotate a well defined cell type. Defaults to -0.2.
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
        max_considered_genes: int = 1000,
        save_output: bool = False,
        input_prefix: str = "aggregated",
        output_prefix: str = "batch_corrected",
        verbose: bool = False) -> list[str]:
    """
    Runs Batch_correction.py on an aggregated h5ad file to correct batch effects.

    Input should be an aggregated h5ad file.

    Outputs a gzip compressed corrected h5ad file with batch corrected representations.
    Output files are named {output_prefix}_{basename}.h5ad.
    
    Annotations added to adata.layers: [dtype(adata.X): "X_log_normalized"] (log1p normalized (total to 1e4) representation of adata.X)
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
        - adata.obs["cancer_state"] (from batch_aggregation.py) (for root cell selection)

    Outputs a gzip compressed corrected h5ad file with CNV annotations.
    Output files are named {output_prefix}_{basename}.h5ad.

    Annotations added to adata.obs: [numpy.float64:"summed_cnvs", numpy.float64:"cnv_score"] (number of CNVs in a cell, cnv score as computed by inferCNVpy.tl.cnv_score)
    Annotations added to adata.obsm: [scipy.sparse._csr.csr_matrix[numpy.float64]:"{corrected_representation}_cnv", numpy.ndarray[numpy.float64]:"{corrected_representation}_gene_values_cnv"] 
    (smoothed and denoised gene expression along genomic locations (basically "cnvness of a gene"), number of copies per gene)

    Paramaeters:
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

    Return:
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
    os.rename(temp_output_path, os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}"))
    temp_output_path = os.path.join(output_temp_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix + "_")}")

    # add output file to output_file_list
    output_file_list.append(temp_output_path)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        print(f"Saving {temp_output_path} to {output_storage_dir}")
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


def infer_pseudotime(
        input_data_file: str, 
        corrected_representation: str = None,
        save_output: bool = False,
        input_prefix: str = "CNV_inferred",
        output_prefix: str = "pseudotime_inferred",
        verbose: bool = False) -> list[str]:
    """ 
    Infer pseudotimefrom a CNV inferred h5ad file. Chooses the cell with least CNVs as root cell.

    Input should be a CNV inferred h5ad file. CNVs are necessary for root cell selection.

    Requires the following annotations to be present:
        - adata.obsm[f"{corrected_representation}_gene_values_cnv"] (from infer_CNV.py)

    Outputs a gzip compressed h5ad with pseudotime annotated for each cell. 
    Output files are named {output_prefix}_{basename}.h5ad.

    Annotations added to adata.obs: adata.obs['dpt_pseudotime']: pandas.Series (dtype float) (pseudotime value for each cell)

    Parameters:
        input_data_file (str): path to aggregated / batch corrected h5ad file.
        corrected_representation (str, optional): name of corrected representation that was used to infer CNVs.
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
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "pseudotime_inference.py"), input_data_file, output_temp_dir, [corrected_representation, verbose])

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

def cluster_and_plot(
        # necessary arguments
        input_data_file: str, 
        modules: list[Literal["projections", "pseudotime_vs_cnv"]],

        # general arguments
        cell_type: str = None, 

        # UMAP / PCA plotting arguments
        obs_annotations: list[str] = None, 
        layers: list[str] = ["X"],
        projection: Literal["UMAP", "PCA"]  = "UMAP",
        marker_file_path: str = None,
        root_cell_idx: int = None,

        # auxiliary arguments
        show: bool = True,
        save_output: bool = False,
        verbose: bool = False) -> None:
    
    """
    Runs plotting modules passed in modules list.

    Current modules include:
        - projections, which will produce UMAP or PCA plots computed from the passed layer
          and colored by the passed obs annotations + leiden clusters. Also computes DEGs for categorical obs annotations + leiden clusters.
          This module uses the following arguments: obs_annotations, layers, projection, marker_file_path, root_cell_idx
        - pseudotime_vs_cnv, which will produce a plot scatterplot of cells, with cnvs on the y axis and pseudotime on the x axis.
          This module uses the following arguments: None

    Input should be an h5ad file with obs annotations for each cell that should be used for clustering.

    Annotations that need to be present:
        - if cell_type is passed, adata.obs["cell_type"] (from cell_type_annotation.py)
        - a column in adata.obs for each passed obs annotation
        - the specified layers need to exist in adata.layers or adata.obsm
        - if modules includes "pseudotime_vs_cnv", adata.obs["dpt_pseudotime"] (from pseudotime_inference.py), adata.obs["cnv_score"], adata.obs["summed_cnvs] (from infer_CNV.py)

    Outputs png files of plots and text files with DEGs if save_output is True. Shows plots if show is True.

    Does not add annotations, as it does not return an h5ad file.

    Parameters:
        input_data_dir (str): path to h5ad file with obs annotations for each cell.
        modules (list[Literal["projections", "pseudotime_vs_cnv"]]): list of modules in plotting scirpt to run.
        
        cell_type (str, optional): Name of the cell type in adata.obs["cell_type"] that should be isolated before plotting.

        obs_annotations (list[str], optional): list of obs annotations to cluster. eg ["cell_type", "pseudotime", "cnv_score", ...].
        layers (str, optional): Anndata layer or obsm keys to use as base for UMAP / PCA projection. Eg. ["X_cnv", "X_scvi_corrected", ...]. Defaults to ["X"], meaning adata.X will be used. 
        projection (literal ["UMAP", "PCA"], optional): projection to use for clustering. Defaults to "UMAP".
        marker_file_path (str, optional): path to file with marker genes for each cell type. DEG analysis will not yield matrix plots if not provided.
        root_cell_idx (int, optional): Index of the root cell of pseudotime. Will be highlighted in pseudotime colored plot.
        
        show (bool, optional): whether to show the plots. Defaults to True.
        save_output (bool, optional): whether to save the plots permanently. Defaults to False.
        verbose (bool, optional): whether to print verbose output. Defaults to False.

    Returns:
        None, but saves files with plots if save_output is True
    """
    # sanity check
    if show == False and save_output == False:
        raise ValueError("At least one of show or save_output must be True")

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have relevant folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "plots"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "plots"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "plots")
    output_temp_dir = os.path.join(TEMP_DIR, "plots")

    # run script and assign path to temporary output file
    print(f"Clustering and plotting: {input_data_file}")
    hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Plotting.py"), input_data_file, output_temp_dir, [modules, cell_type, obs_annotations, layers, projection, marker_file_path, root_cell_idx, show, verbose])

    # naming happens in subprocess (relies on knowing which obs column was used)

    # if specified, permanently store a copy of the temporary output file
    if save_output == True:
        for file in os.listdir(output_temp_dir):
            shutil.copy(os.path.join(output_temp_dir, file), os.path.join(output_storage_dir, file))

    # outputs will not be used programatically so no need to return list with output file paths
    # this is why we also don't use input and output prefixes
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

        #mode = choose_pipeline_mode(RAW_DATA_DIRS[0])
        # temp_output_files = preprocess_data(RAW_DATA_DIRS[0], mode, use_ensembl_ids=use_ensebml_ids, save_output=False, verbose=True, filtering_params=FilteringParameters(min_n_cells_percentage=0))
        # annotate_cell_types(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed"), use_ensebml_ids, r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\marker_genes.json", r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\negative_markers.json", model="cellassign", verbose=True, save_output=True)
        output_path_list = aggregate_batches(os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated"), save_output=False, verbose=True)
        correct_batch_effects(output_path_list[0], save_output=True, verbose=True, max_considered_genes="all")
        # infer_CNVs(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\aggregated\aggregated_PDAC.h5ad", r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\gencode.v49.annotation.gtf.gz", save_output=True, input_prefix="aggregated", verbose=True, cell_type="ductal_cell")
        # infer_pseudotime(os.path.join(OUTPUT_STORAGE_DIR, "CNV", "CNV_inferred_PDAC.h5ad"), verbose=True, corrected_representation=None, save_output=True)
        # cluster_and_plot(os.path.join(OUTPUT_STORAGE_DIR, "pseudotime", "pseudotime_inferred_PDAC.h5ad"), ["projections", "pseudotime_vs_cnv"], layers=["X", "X_cnv"], root_cell_idx=14, marker_file_path=r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\marker_genes.json", obs_annotations=["cnv_score", "dpt_pseudotime", "cell_type", "summed_cnvs", "cancer_state", "batch"], cell_type="ductal_cell", projection="UMAP", show=True, save_output=False, verbose=True)

        purge_tempfiles()
        sys.exit(0)
    except Exception:
        purge_tempfiles()
        raise # re-raise the exception to see the traceback and error message


