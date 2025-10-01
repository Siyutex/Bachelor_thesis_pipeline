# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients
# change RAW_DATA_DIR and outcome_storage for each run to reflect new input and what outcomes should be permanenlty saved

import subprocess # needed for running other scripts
import os # needed for file and directory operations
from enum import Enum
import json # needed for forwarding python objects to subprocesses
import shutil # needed for file storage operations
import tempfile # needed for temporary file operations
import sys # needed to exit the program
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
                os.path.join(SCRIPT_DIR, "..", "..", "Data","pretraining", "cancerSCEM", "colon_cancer_cancerous"),    
                os.path.join(SCRIPT_DIR, "..", "..", "Data","OG_data","NCBI","PDAC_cancerous"),
                ]
OUTPUT_STORAGE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data", "output_storage")  # directory for optional permanent storage of indermediate subprocess outputs
TEMP_DIR = os.path.join(tempfile.gettempdir(),"python") # directory for storage of temporary pipeline files


# variable to determine what intermediate files should be saved permanently, one key per script
OUTCOME_STORAGE = {
    "Preprocessing.py": True,
    "Cell_type_annotation.py": True,
    "Clustering.py": True,
    "Batch_correction.py": True,
    "GRN_inference.py": True,

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



def preprocess_data(pipeline_mode: pipeline_mode, raw_data_dir: str):
    """Loops through the RAW_DATA_DIR and runs Preprocessing.py on each file / folder based on the chosen pipeline mode.
    Saves preprocessed files to temp/preprocessed if outcome_storage for Preprocessing.py is True."""

    # check if outcome storage directory has preprocessed folder, if not create it
    if not os.path.exists(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed")) and OUTCOME_STORAGE["Preprocessing.py"] == True:  # check if temp directory has preprocessed folder and outcome storage is True
        os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed"))  # create preprocessed folder if it does not exist 
    
    # output path for preprocessed files with a prefix in a subdirectory of OUTCOME_STORAGE_DIR
    if OUTCOME_STORAGE["Preprocessing.py"] == True:
        output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "preprocessed")

    # check if preprocessed folder exists in TEMP_DIR, if not create it
    os.makedirs(os.path.join(TEMP_DIR, "preprocessed"), exist_ok=True)

    # assign output path variable to be equal to TEMP_DIR/preprocessed
    output_temp_dir = os.path.join(TEMP_DIR, "preprocessed")

    # assign datatype variable based on chosen pipeline mode
    datatype = str(pipeline_mode.name)

    # run script on RAW_DATA_DIR based on chosen pipeline mode
    # define iterator for file naming + 1 for each file in element in raw_data_dir
    i = 0
    if datatype == pipeline_mode.MTX_TSVs_in_subfolders.name:
        
        # iterate over folders in raw data directory containing two tsv files and one mtx file each
        for folder in os.listdir(raw_data_dir):
            print("Currently preprocessing: " +  folder)
            temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, folder), output_temp_dir, [datatype])
            
            # rename file at temp_output_path to "preprocessed_{raw_data_dir}_i.h5ad" and adjust path
            os.rename(temp_output_path, os.path.join(output_temp_dir, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
            temp_output_path = os.path.join(output_temp_dir, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad")

            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

            i += 1

    elif datatype == pipeline_mode.compressed_MTX_TSVs_in_subfolders.name:

        # iterate over folder in raw data directory, then forward compressed files in them
        for folder in os.listdir(raw_data_dir):
            print("Currently preprocessing: " +  folder)
            temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, folder), output_temp_dir, [datatype])

            # rename file at temp_output_path to "preprocessed_{raw_data_dir}_i.h5ad" and adjust path
            os.rename(temp_output_path, os.path.join(output_temp_dir, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
            temp_output_path = os.path.join(output_temp_dir, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad")

            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

            i += 1

    elif datatype == pipeline_mode.dot_matrix_files.name:

        # directly forward .matrix files in RAW_DATA_DIR
        for file in os.listdir(raw_data_dir):
            print("Currently preprocessing: " +  file)
            temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, file), output_temp_dir, [datatype])

            # rename file at temp_output_path to "preprocessed_{raw_data_dir}_i.h5ad" and adjust path
            os.rename(temp_output_path, os.path.join(output_temp_dir, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
            temp_output_path = os.path.join(output_temp_dir, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad")

            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

            i += 1


def annotate_cell_types(input_data_dir: str):
    """ 
    Run Cell_type_annotation.py on a given directory of preprocessed h5ad files.
    Creates necessary directories if not present. Saves output permanenlty if specified.
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have cell_type_annotated folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "cell_type_annotated"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "cell_type_annotated")
    output_temp_dir = os.path.join(TEMP_DIR, "cell_type_annotated")

    # run script on each file in input_data_dir
    for file in os.listdir(input_data_dir):
        print("Annotating cell types for: " + file)
        temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "Cell_type_annotation.py"), os.path.join(input_data_dir, file), output_temp_dir)

        # if specified, permanently store a copy of the temporary output file
        if OUTCOME_STORAGE["Cell_type_annotation.py"] == True:
            shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))


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


def infer_GRN(input_data_dir: str, verbose: bool = False):
    """ 
    Infer a GRN (csv, regulators + importance) from a given aggregated / batch corrected h5ad file.
    """

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "GRN"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "GRN"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "GRN")
    output_temp_dir = os.path.join(TEMP_DIR, "GRN")

    # run script and assign path to temporary output file
    print(f"Inferring GRN from {input_data_dir}")
    temp_output_path = hf.execute_subprocess(os.path.join(SCRIPT_DIR, "GRN_inference.py"), input_data_dir, output_temp_dir, [verbose])

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["GRN_inference.py"] == True:
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

    try:
        '''for raw_data_dir in RAW_DATA_DIRS:
            mode = choose_pipeline_mode(raw_data_dir)
            preprocess_data(mode, raw_data_dir)'''
        for file in os.listdir(os.path.join(OUTPUT_STORAGE_DIR,"batch_corrected")):
            infer_GRN(os.path.join(OUTPUT_STORAGE_DIR,"batch_corrected",file), verbose=True)
        purge_tempfiles()
        sys.exit(0) # don't want to loop, while is just to be able to break out of it with a signal
    except Exception:
        purge_tempfiles()
        raise # re-raise the exception to see the traceback and error message
    




