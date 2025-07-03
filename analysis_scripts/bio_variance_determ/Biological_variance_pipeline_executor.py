# This script runs the pipeline to preprocess mtx + tsv files from the 10x genomics pipeline and isolate epithelial cells from them and determine biological variance between patients

import subprocess # needed for running other scripts
import os # needed for file and directory operations

#check if the required directories exist, if not create them
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data")):  # check if data directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data"))  # create data directory if it does not exist
    print("Created Data directory. Please add subdirectories for each sample containing the mtx and tsv outputs of the 10x genomics pipeline")  # inform user to add subdirectories with required files
if not os.path.exists(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "temp")):  # check if temp directory exists
    os.makedirs(os.path.join(os.path.dirname(__file__), "..", "..", "Data", "temp"))  # create temp directory if it does not exist
    print("Created temp directory. Preprocessed files will be saved here.")  # inform user that preprocessed files will be saved here


# declare script directory and mtx directory as global constants
script_dir = os.path.dirname(__file__)  # directory where this script is located
raw_data_dir = os.path.join(script_dir, "..", "..", "Data")  # directory where mtx files are located
temp_dir = os.path.join(script_dir, "..", "..", "Data", "temp")  # directory for temporary files

# functions
def preprocess_data():
    """Loops through folders in data directory and runs Preprocessing.py on each. 
    The folders must each contain "matrix.mtx", "genes.tsv" and "barcodes.tsv".
    The output is saved in the temp directory."""
    
    # check if the data directory contains subdirectories with mtx files
    if not any(any(file.endswith(".mtx") for file in os.listdir(os.path.join(raw_data_dir, folder))) for folder in os.listdir(raw_data_dir)): # check if any subdirectory contains mtx files
        print("Preprocessing: No subdirectories with mtx files found in the data directory. Please add subdirectories with mtx and tsv files from 10x genomics pipeline.")
        return # exit the function if no subdirectories with mtx files are found
    
    # check if temp directory has preprocessed folder, if not create it
    if not os.path.exists(os.path.join(temp_dir, "preprocessed")):  # check if temp directory has preprocessed folder
        os.makedirs(os.path.join(temp_dir, "preprocessed"))  # create preprocessed folder if it does not exist
    
    # iterate over folders in raw data directory containing two tsv files and one mtx file each
    for folder in os.listdir(raw_data_dir): 
        # skip temp folder, folders that do not contain mtx files and check if the folder is a directory
        if os.path.join(raw_data_dir, folder) != os.path.join(raw_data_dir, "temp"):
            print("Currently preprocessing: " +  folder)
            output_path = os.path.join(script_dir, "..", "..", "Data", "temp", "preprocessed", f"preprocessed_{folder}")  # output path for preprocessed files with a prefix in a subdirectory of temp directory
            subprocess.run(["python", os.path.join(script_dir, "Preprocessing.py"), os.path.join(raw_data_dir, folder), output_path], check=True) # check=True ensures that an error in the subprocess will raise an exception
        else:
            print(f"Preprocessing: Skipping folder {folder} as it does not contain mtx files or is the temp directory.")



def isolate_epithelial_cells():
    """Loops through the preprocessed h5ad files in temp/preprocessed and runs Epithelial_cell_isolation.py on each."""

    for file in os.listdir(os.path.join(temp_dir, "preprocessed")):  # iterate over files in the temp directory where preprocessed files are stored
        if file.endswith(".h5ad"): # check if the file is an h5ad file
            (print("Currently isolating epithelial cells from: " + file))
            file_path = os.path.join(temp_dir, "preprocessed", file)
            output_path = os.path.join(temp_dir, f"epithelial_isolated_{file}") # save output in the temp directory with another new prefix
            subprocess.run(["python", os.path.join(script_dir, "Epithelial_cell_isolation.py"), file_path, output_path], check=True) # check=True ensures that an error in the subprocess will raise an exception
        else:
            print(f"Epithelial_isolation: Skipping file {file} as it is not a .h5ad file.")

if __name__ == "__main__": # ensures this code runs only when this script is executed directly, not when imported
    preprocess_data()
    isolate_epithelial_cells()