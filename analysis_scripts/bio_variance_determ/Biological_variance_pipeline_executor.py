# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients
# change RAW_DATA_DIR and outcome_storage for each run to reflect new input and what outcomes should be permanenlty saved

import subprocess # needed for running other scripts
import os # needed for file and directory operations
from enum import Enum
import json # needed for forwarding python objects to subprocesses
import shutil # needed for file storage operations
import tempfile # needed for temporary file operations
import sys # needed to exit the program

from rich.traceback import install # needed for pretty printing of tracebacks
install(show_locals=True, suppress=[]) # install rich traceback handler


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
                os.path.join(SCRIPT_DIR, "..", "..", "Data","pretraining","cancerSCEM","clear_cell_renal_carcinoma_cancerous"),
                os.path.join(SCRIPT_DIR, "..", "..", "Data","pretraining","cancerSCEM","clear_cell_renal_carcinoma_non_cancerous"),
                ]
OUTPUT_STORAGE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data", "output_storage")  # directory for optional permanent storage of indermediate subprocess outputs
TEMP_DIR = os.path.join(tempfile.gettempdir(),"python") # directory for storage of temporary pipeline files


# variable to determine what intermediate files should be saved permanently, one key per script
OUTCOME_STORAGE = {
    "Preprocessing.py": False,
    "Batch_correction.py": False,
    "Epithelial_cell_isolation.py": False,
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


# temporary comment: arguments will be: subprocess path, rawdata path, output path (not needed, bcs output is handled in executor), datatype, list of python objects to be forwarded, outcome storage dictionary (not needed, bcs handled per function that runs a subprocess)
import subprocess
import json
import os

def execute_subprocess(subprocess_path: str, inputadata_path: str, output_dir_path: str,
                       inputdatatype: str = None, python_objects: list = None) -> str:
    """Run a Python subprocess, forward arguments and optional JSON objects, return output file path."""

    args = ["python", "-u", subprocess_path, inputadata_path, output_dir_path, inputdatatype or "", str(len(python_objects or []))] # u means unbuffered mode -> directly print when print statement

    # Prepare stdin payload if needed
    stdin_data = None
    if python_objects:
        stdin_data = "\n".join(json.dumps(obj) for obj in python_objects) + "\n"

    # Run subprocess, capture stdout + stderr
    proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Parse every line of stdout
    output_file_path = None
    for line in proc.stdout:
        if line.strip().startswith("Output: "):
            output_file_path = line.split("Output:", 1)[1].strip()
            print(f"execute_subprocess: received output: {output_file_path}")
        else:
            print(f"Subprocess {os.path.basename(subprocess_path)}: {line}")

    # Wait for subprocess to finish and collect stderr
    stdout, stderr = proc.communicate(input=stdin_data)

    # if there was an error, raise RntimeError and print stderr
    if proc.returncode != 0:
        raise RuntimeError(f"Subprocess {os.path.basename(subprocess_path)} failed with exit code {proc.returncode}:\n{stderr}")

    # if output file path was not printed to stdout by the subprocess, raise RuntimeError
    if output_file_path is None:
        raise RuntimeError(f"Output file path not found in subprocess {os.path.basename(subprocess_path)} stdout")

    return output_file_path




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
        output_storage_path = os.path.join(OUTPUT_STORAGE_DIR, "preprocessed")

    # check if preprocessed folder exists in TEMP_DIR, if not create it
    os.makedirs(os.path.join(TEMP_DIR, "preprocessed"), exist_ok=True)

    # assign output path variable to be equal to TEMP_DIR/preprocessed
    output_path = os.path.join(TEMP_DIR, "preprocessed")

    # assign datatype variable based on chosen pipeline mode
    datatype = pipeline_mode.name

    # run script on RAW_DATA_DIR based on chosen pipeline mode
    # define iterator for file naming + 1 for each file in element in raw_data_dir
    i = 0
    if datatype == pipeline_mode.MTX_TSVs_in_subfolders.name:
        
        # iterate over folders in raw data directory containing two tsv files and one mtx file each
        for folder in os.listdir(raw_data_dir):
            print("Currently preprocessing: " +  folder)
            temp_output_path = execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, folder), output_path, datatype)
            
            # rename file at temp_output_path to "preprocessed_{raw_data_dir}_i.h5ad" and adjust path
            os.rename(temp_output_path, os.path.join(output_path, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
            temp_output_path = os.path.join(output_path, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad")

            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_storage_path, os.path.basename(temp_output_path)))

            i += 1

    elif datatype == pipeline_mode.compressed_MTX_TSVs_in_subfolders.name:

        # iterate over folder in raw data directory, then forward compressed files in them
        for folder in os.listdir(raw_data_dir):
            print("Currently preprocessing: " +  folder)
            temp_output_path = execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, folder), output_path, datatype)

            # rename file at temp_output_path to "preprocessed_{raw_data_dir}_i.h5ad" and adjust path
            os.rename(temp_output_path, os.path.join(output_path, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
            temp_output_path = os.path.join(output_path, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad")

            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_storage_path, os.path.basename(temp_output_path)))

            i += 1

    elif datatype == pipeline_mode.dot_matrix_files.name:

        # directly forward .matrix files in RAW_DATA_DIR
        for file in os.listdir(raw_data_dir):
            print("Currently preprocessing: " +  file)
            temp_output_path = execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(raw_data_dir, file), output_path, datatype)

            # rename file at temp_output_path to "preprocessed_{raw_data_dir}_i.h5ad" and adjust path
            os.rename(temp_output_path, os.path.join(output_path, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad"))
            temp_output_path = os.path.join(output_path, f"preprocessed_{os.path.basename(raw_data_dir)}_{i}.h5ad")

            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_storage_path, os.path.basename(temp_output_path)))

            i += 1


def correct_batch_effects(input_data_dir: str):
    """Run Batch_correction.py on a given directory of preprocessed h5ad files."""

    #check if OUTCOME_STORAGE_DIR and TEMP_DIR have batch_corrected folder, if not create it
    os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected"), exist_ok=True)
    os.makedirs(os.path.join(TEMP_DIR, "batch_corrected"), exist_ok=True)

    # assign directories for temporary and permanent storage
    output_storage_dir = os.path.join(OUTPUT_STORAGE_DIR, "batch_corrected")
    output_temp_dir = os.path.join(TEMP_DIR, "batch_corrected")

    # run script and assign path to temporary output file
    print(f"Correcting batch effects in {input_data_dir}")
    temp_output_path = execute_subprocess(os.path.join(SCRIPT_DIR, "Batch_correction.py"), input_data_dir, output_temp_dir)

    # if specified, permanently store a copy of the temporary output file
    if OUTCOME_STORAGE["Batch_correction.py"] == True:
        shutil.copy(temp_output_path, os.path.join(output_storage_dir, os.path.basename(temp_output_path)))

def isolate_epithelial_cells():
    """Loops through the preprocessed h5ad files in temp/preprocessed and runs Epithelial_cell_isolation.py on each."""

    # check if temp directory has epithelial_isolated folder, if not create it
    if not os.path.exists(os.path.join(OUTPUT_STORAGE_DIR, "epithelial_isolated")):  # check if temp directory has epithelial_isolated folder
        os.makedirs(os.path.join(OUTPUT_STORAGE_DIR, "epithelial_isolated"))  # create epithelial_isolated folder if it does not exist

    for file in os.listdir(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed")):  # iterate over files in the temp directory where preprocessed files are stored
        if file.endswith(".h5ad"): # check if the file is an h5ad file
            (print("Currently isolating epithelial cells from: " + file))
            file_path = os.path.join(OUTPUT_STORAGE_DIR, "preprocessed", file)
            output_path = os.path.join(OUTPUT_STORAGE_DIR, "epithelial_isolated", f"epithelial_isolated_{file}") # save output in the temp directory with another new prefix
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Epithelial_cell_isolation.py"), file_path, output_path], check=True) # check=True ensures that an error in the subprocess will raise an exception
        else:
            print(f"Epithelial_isolation: Skipping file {file} as it is not a .h5ad file.")


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
        correct_batch_effects(os.path.join(OUTPUT_STORAGE_DIR, "preprocessed"))
        purge_tempfiles()
        sys.exit(0) # don't want to loop, while is just to be able to break out of it with a signal
    except Exception:
        purge_tempfiles()
        raise # re-raise the exception to see the traceback and error message
    




