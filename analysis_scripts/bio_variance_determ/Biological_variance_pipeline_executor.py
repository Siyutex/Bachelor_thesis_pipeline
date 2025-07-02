# This script runs the pipeline to preprocess mtx files and isolate epithelial cells from them and determine biological variance between patients

import subprocess # needed for running other scripts
import os # needed for file and directory operations

# declare script directory and mtx directory as global constants
script_dir = os.path.dirname(__file__)  # directory where this script is located
raw_data_dir = os.path.join(script_dir, "..", "..", "Data", "PDAC_vs_P_Epi_scRNAseq")  # directory where mtx files are located
temp_dir = os.path.join(script_dir, "..", "..", "Data", "temp")  # directory for temporary files

# functions
def preprocess_data():
    """Loops through folders in data directory and runs Preprocessing.py on each. 
    The folders must each contain "matrix.mtx", "genes.tsv" and "barcodes.tsv".
    The output is saved in the temp directory."""
    
    for folder in os.listdir(raw_data_dir): # iterate over folders in raw data directory containing two tsv files and one mtx file each
        # skip temp folder, folders that do not contain mtx files and check if the folder is a directory
        if os.path.join(raw_data_dir, folder) != os.path.join(raw_data_dir, "temp") and any(file.endswith(".mtx") for file in os.listdir(os.path.join(raw_data_dir, folder))):
            print("Currently preprocessing: " + os.path.join(raw_data_dir, folder))
            output_path = os.path.join(script_dir, "..", "..", "Data", "PDAC_vs_P_Epi_scRNAseq", "temp", f"preprocessed_{folder}")  # output path for preprocessed files with a prefix in the temp directory, forwarded to Preprocessing.py
            subprocess.run(["python", os.path.join(script_dir, "Preprocessing.py"), os.path.join(raw_data_dir, folder), output_path], check=True) # check=True ensures that an error in the subprocess will raise an exception
        else:
            print(f"Preprocessing: Skipping folder {folder} as it does not contain mtx files or is the temp directory.")



def isolate_epithelial_cells():
    """Loops through the preprocessed h5ad files in temp and runs Epithelial_cell_isolation.py on each."""

    for file in os.listdir(temp_dir):  # iterate over files in the temp directory where preprocessed files are stored
        if file.endswith(".mtx"):
            file_path = os.path.join(temp_dir, file)
            output_path = os.path.join(temp_dir, f"epithelial_isolated_{file}") # save output in the temp directory with another new prefix
            subprocess.run(["python", os.path.join(script_dir, "Epithelial_cell_isolation.py"), file_path, output_path], check=True) # check=True ensures that an error in the subprocess will raise an exception
        elif file.endswith(".mtx") != True: 
            print(f"Epithelial_isolation: Skipping file {file} as it is not a .mtx file.")

if __name__ == "__main__": # ensures this code runs only when this script is executed directly, not when imported
    preprocess_data()