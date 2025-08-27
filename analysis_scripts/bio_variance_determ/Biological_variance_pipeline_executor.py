# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients
# change RAW_DATA_DIR and outcome_storage for each run to reflect new input and what outcomes should be permanenlty saved

import subprocess # needed for running other scripts
import os # needed for file and directory operations
from enum import Enum
import json # needed for forwarding python objects to subprocesses
import shutil # needed for file storage operations
import tempfile # needed for temporary file operations
import time

#check if the required directories exist, if not create them
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data")):  # check if data directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data"))  # create data directory if it does not exist
    print("Created Data directory. Please add subdirectories for each sample containing the mtx and tsv outputs of the 10x genomics pipeline")  # inform user to add subdirectories with required files
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "temp")):  # check if temp directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "temp"))  # create temp directory if it does not exist
    print("Created temp directory. Preprocessed files will be saved here.")  # inform user that preprocessed files will be saved here


# declare script directory and mtx directory as global constants
SCRIPT_DIR = os.path.dirname(__file__)  # directory where this script is located
RAW_DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data","OG_PDAC_and_Pancreas")  # directory where files / folder with files are located (10x genomics, GDC or cancerSCEM format)
TEMP_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data", "temp")  # directory for temporary files

# variable to determine what intermediate files should be saved permanently, one key per script
OUTCOME_STORAGE = {
    "Preprocessing.py": True,
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

def choose_pipeline_mode():
    """Inspects the structure of the provided RAW_DATA_DIR and returns the appropraote
    pipeline_mode enum value."""
    
    mode = pipeline_mode.NO_MODE_CHOSEN

    # check if RAW_DATA_DIR has subdirectories (to avoid errors), if so check their structure
    if any(os.path.isdir(os.path.join(RAW_DATA_DIR, folder)) for folder in os.listdir(RAW_DATA_DIR)):

        if any(any(file.endswith(".mtx") for file in os.listdir(os.path.join(RAW_DATA_DIR,folder))) for folder in os.listdir(RAW_DATA_DIR)):
            # if there is an mtx in any subfolder of RAW_DATA_DIR, we assume that it is full of folders with mtx + 2 tsvs (10x genomics format)
            mode = pipeline_mode.MTX_TSVs_in_subfolders

        if any(any(file.endswith("matrix.tar.gz") for file in os.listdir(os.path.join(RAW_DATA_DIR,folder))) for folder in os.listdir(RAW_DATA_DIR)):
            # if there is such a structure, but with .gz compression (GDC format)
            mode = pipeline_mode.compressed_MTX_TSVs_in_subfolders
    elif any(file.endswith("counts.matrix.tsv.gz") for file in os.listdir(RAW_DATA_DIR)):
        # if the RAW_DATA_DIR directly contains .matrix files (cancerSCEM format)
        mode = pipeline_mode.dot_matrix_files



    if mode == pipeline_mode.NO_MODE_CHOSEN:
        raise ValueError("Could not determine pipeline mode. Please check the structure of the provided RAW_DATA_DIR. It should either contain subdirectories with mtx and tsv files (10x genomics format), subdirectories with compressed mtx and tsv files (GDC format) or directly .matrix files (cancerSCEM format).")

    print(f"Pipeline mode: {mode.name}")
    return mode


# temporary comment: arguments will be: subprocess path, rawdata path, output path (not needed, bcs output is handled in executor), datatype, list of python objects to be forwarded, outcome storage dictionary (not needed, bcs handled per function that runs a subprocess)
def execute_subprocess(subprocess_path: str, inputadata_path: str, inputdatatype: str=None, python_objects: list = None):
    """Runs a python subprocess, forwards python objects to it, returns a path to a temporary file with the output,
    optionally saves output to file. Only subprocess_path and inputadata_path are required arguments. 
    Inputadatatype is only needed for Preprocessing.py."""

    if python_objects is not None:
        num_objects = str(len(python_objects))  # number of python objects to be forwarded, needed in subprocess to know how many to expect
    else:
        num_objects = "0"  # if no python objects are to be forwarded, set to 0 (needs to be a string to be passed as command line argument)
    
    # intialize the subprocess and pass command line arguments (this starts the subprocess, but you can still interact with it (not possible with subprocess.run))
    # subprocess.run and subprocess.Popen take args = [] as first argument, args[0] is the executable (the python itnerpreter), and the rest are arguments to the interpreter (so args[1], in this case is a filepath to a python executable, so it will be executed by the itnerpreter)
    proc = subprocess.Popen(["python", subprocess_path, inputadata_path, inputdatatype, num_objects], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # forward python objects to subprocess via stdin (need to be read in the subprocess)
    if python_objects is not None:
        for python_object in python_objects:
            json_string = json.dumps(python_object)  # convert python object to json string
            proc.stdin.write(json_string + "\n")  # write json string to stdin of subprocess, followed by newline
            proc.stdin.flush()  # flush the stdin buffer to ensure the data is sent

    # send back the path to a temporary file with the output
    # temp files should be tempfile.NamedTemporaryFile(delete=False) so they still exist after the subprocess ends
    for line in proc.stdout:
        # we expect the paths that are returned from the subprocesses to be structured like "Output: " + [path]
        if line.strip().startswith("Output: "):
            output_path = line.split("Output:")[1]  # read the output from stdout of subprocess (should be a path to a temporary file)
            print(f"execute_subrpocess: recieved output: {output_path}")
            return output_path.strip()  # return the path to the temporary file with the output


def purge_tempfiles():
    # get the system temp directory (e.g. /tmp on Linux, AppData\Local\Temp on Windows)
    os.makedirs(os.path.join(tempfile.gettempdir(),"python"), exist_ok=True)
    pipeline_temp_dir = os.path.join(tempfile.gettempdir(),"python")

    # careful: this deletes *everything* inside that dir
    for filename in os.listdir(pipeline_temp_dir):
        file_path = os.path.join(pipeline_temp_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)       # delete file or symlink
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)   # delete folder and contents
        except Exception as e:
            print(f"Failed to delete {file_path}: {e}")

def preprocess_data(pipeline_mode: pipeline_mode):
    """Loops through the RAW_DATA_DIR and runs Preprocessing.py on each file / folder based on the chosen pipeline mode.
    Saves preprocessed files to temp/preprocessed if outcome_storage for Preprocessing.py is True."""

    # check if temp directory has preprocessed folder, if not create it
    if not os.path.exists(os.path.join(TEMP_DIR, "preprocessed")) and OUTCOME_STORAGE["Preprocessing.py"] == True:  # check if temp directory has preprocessed folder and outcome storage is True
        os.makedirs(os.path.join(TEMP_DIR, "preprocessed"))  # create preprocessed folder if it does not exist 
    
    # output path for preprocessed files with a prefix in a subdirectory of TEMP_DIR
    if OUTCOME_STORAGE["Preprocessing.py"] == True:
        output_path = os.path.join(TEMP_DIR, "preprocessed")

    # assign datatype variable based on chosen pipeline mode
    datatype = pipeline_mode.name
    print(f"Chosen pipeline mode: {datatype}")
    print(pipeline_mode.MTX_TSVs_in_subfolders.name)

    # run script on RAW_DATA_DIR based on chosen pipeline mode
    if datatype == pipeline_mode.MTX_TSVs_in_subfolders.name:
        
        # iterate over folders in raw data directory containing two tsv files and one mtx file each
        for folder in os.listdir(RAW_DATA_DIR):
             
            print("Currently preprocessing: " +  folder)

            temp_output_path = execute_subprocess(os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(RAW_DATA_DIR, folder), datatype)
            
            # if specified, permanently store a copy of the temporary output file
            if OUTCOME_STORAGE["Preprocessing.py"] == True:
                shutil.copy(temp_output_path, os.path.join(output_path, f"preprocessed_{folder}.h5ad"))
    
    elif datatype == pipeline_mode.compressed_MTX_TSVs_in_subfolders.name:
        # iterate over folder in raw data directory, then forward compressed files in them
        for folder in os.listdir(RAW_DATA_DIR):
            print("Currently preprocessing: " +  folder)
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(RAW_DATA_DIR, folder), os.path.join(output_path, folder), datatype, OUTCOME_STORAGE], check=True) # check=True ensures that an error in the subprocess will raise an exception
    elif datatype == pipeline_mode.dot_matrix_files.name:
        # directly forward .matrix files in RAW_DATA_DIR
        for file in os.listdir(RAW_DATA_DIR):
            print("Currently preprocessing: " +  file)
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(RAW_DATA_DIR, file), os.path.join(output_path, file), datatype, OUTCOME_STORAGE], check=True) # check=True ensures that an error in the subprocess will raise an exception

def isolate_epithelial_cells():
    """Loops through the preprocessed h5ad files in temp/preprocessed and runs Epithelial_cell_isolation.py on each."""

    # check if temp directory has epithelial_isolated folder, if not create it
    if not os.path.exists(os.path.join(TEMP_DIR, "epithelial_isolated")):  # check if temp directory has epithelial_isolated folder
        os.makedirs(os.path.join(TEMP_DIR, "epithelial_isolated"))  # create epithelial_isolated folder if it does not exist

    for file in os.listdir(os.path.join(TEMP_DIR, "preprocessed")):  # iterate over files in the temp directory where preprocessed files are stored
        if file.endswith(".h5ad"): # check if the file is an h5ad file
            (print("Currently isolating epithelial cells from: " + file))
            file_path = os.path.join(TEMP_DIR, "preprocessed", file)
            output_path = os.path.join(TEMP_DIR, "epithelial_isolated", f"epithelial_isolated_{file}") # save output in the temp directory with another new prefix
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Epithelial_cell_isolation.py"), file_path, output_path], check=True) # check=True ensures that an error in the subprocess will raise an exception
        else:
            print(f"Epithelial_isolation: Skipping file {file} as it is not a .h5ad file.")



def compute_variance():
    """Loops through the epithelial isolated h5ad files in temp/epithelial_isolated and runs Variance.py on each."""

    ADJ_paths = []  # list to store paths to tumor adjacent tissue files
    PDAC_paths = []  # list to store paths to tumor tissue files

    for file in os.listdir(os.path.join(TEMP_DIR, "epithelial_isolated")):  # iterate over files in the temp directory where epithelial isolated files are stored
        if file.endswith(".h5ad") and "ADJ" in file:  # check if the file is a h5ad file and tumor adjacent tissue
            ADJ_paths.append(os.path.join(TEMP_DIR, "epithelial_isolated", file))  # add the file path to the list
        elif file.endswith(".h5ad") and "PDAC" in file:  # check if the file is a h5ad file and tumor tissue
            PDAC_paths.append(os.path.join(TEMP_DIR, "epithelial_isolated", file))
        else:
            print(f"Variance: Skipping file {file} as it is not a .h5ad.gz file.")

    ADJ_paths = ",".join(ADJ_paths)  # join the list of paths into a single string separated by commas (so it can be passed as a command line argument)
    PDAC_paths = ",".join(PDAC_paths)

    subprocess.run(["python", os.path.join(SCRIPT_DIR, "Variance.py"), ADJ_paths, PDAC_paths], check=True)
            



if __name__ == "__main__": # ensures this code runs only when this script is executed directly, not when imported
    mode = choose_pipeline_mode()
    preprocess_data(mode)
    purge_tempfiles()