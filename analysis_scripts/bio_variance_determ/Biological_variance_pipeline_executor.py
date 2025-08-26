# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients

import subprocess # needed for running other scripts
import os # needed for file and directory operations
from enum import Enum

#check if the required directories exist, if not create them
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data")):  # check if data directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data"))  # create data directory if it does not exist
    print("Created Data directory. Please add subdirectories for each sample containing the mtx and tsv outputs of the 10x genomics pipeline")  # inform user to add subdirectories with required files
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "temp")):  # check if temp directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "temp"))  # create temp directory if it does not exist
    print("Created temp directory. Preprocessed files will be saved here.")  # inform user that preprocessed files will be saved here


# declare script directory and mtx directory as global constants
SCRIPT_DIR = os.path.dirname(__file__)  # directory where this script is located
RAW_DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data","pretraining","GDC","Gliomas")  # directory where files / folder with files are located (10x genomics, GDC or cancerSCEM format)
TEMP_DIR = os.path.join(SCRIPT_DIR, "..", "..", "Data", "temp")  # directory for temporary files


# classes
class pipeline_mode(Enum):
    MTX_TSVs_in_subfolders = 1  
    compressed_MTX_TSVs = 2
    dot_matrix_files = 3
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
            mode = pipeline_mode.compressed_MTX_TSVs

    if any(file.endswith("counts.matrix.tsv.gz") for file in os.listdir(RAW_DATA_DIR)):
        # if the RAW_DATA_DIR directly contains .matrix files (cancerSCEM format)
        mode = pipeline_mode.dot_matrix_files

    print(f"Pipeline mode: {mode.name}")

    return mode


def preprocess_data(pipeline_mode: pipeline_mode):
    """Loops through folders in data directory and runs Preprocessing.py on each. 
    The folders must each contain "matrix.mtx", "genes.tsv" and "barcodes.tsv".
    The output is saved in the temp directory."""
        
    # check if temp directory has preprocessed folder, if not create it
    if not os.path.exists(os.path.join(TEMP_DIR, "preprocessed")):  # check if temp directory has preprocessed folder
        os.makedirs(os.path.join(TEMP_DIR, "preprocessed"))  # create preprocessed folder if it does not exist
    
    # assign datatype variable based on chosen pipeline mode
    datatype = pipeline_mode.name
    
    # output path for preprocessed files with a prefix in a subdirectory of temp directory
    output_path = os.path.join(TEMP_DIR, "preprocessed", f"preprocessed_{folder}")

    # run script on RAW_DATA_DIR based on chosen pipeline mode
    if pipeline_mode == pipeline_mode.MTX_TSVs_in_subfolders:
        # iterate over folders in raw data directory containing two tsv files and one mtx file each
        for folder in os.listdir(RAW_DATA_DIR): 
            print("Currently preprocessing: " +  folder)
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(RAW_DATA_DIR, folder), output_path, datatype], check=True) # check=True ensures that an error in the subprocess will raise an exception
    elif pipeline_mode == pipeline_mode.compressed_MTX_TSVs:
        # iterate over folder in raw data directory, then forward compressed files in them
        for folder in os.listdir(RAW_DATA_DIR):
            print("Currently preprocessing: " +  folder)
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(RAW_DATA_DIR, folder), output_path, datatype], check=True) # check=True ensures that an error in the subprocess will raise an exception
    elif pipeline_mode == pipeline_mode.dot_matrix_files:
        # directly forward .matrix files in RAW_DATA_DIR
        for file in os.listdir(RAW_DATA_DIR):
            print("Currently preprocessing: " +  file)
            subprocess.run(["python", os.path.join(SCRIPT_DIR, "Preprocessing.py"), os.path.join(RAW_DATA_DIR, file), output_path, datatype], check=True) # check=True ensures that an error in the subprocess will raise an exception


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