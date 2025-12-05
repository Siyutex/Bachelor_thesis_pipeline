# this script compares which cells (by labels) are assigned to the transition
# state across runs (issue: labelling might be different)

import scanpy as sc
import helper_functions as hf
import os
from collections import Counter
from scipy import sparse
import numpy as np
from collections import defaultdict
import hashlib


def load_obsname_dict(file_list: list[str]):
    """
    takes a list of filepaths to h5ad files
    return a dictionary, key = filepath, values list(adata.obs_names)
    """
    # returns a dict, key = origin file, value = obs_names
    obsname_dict = {}
    for file in file_list:
        adata = sc.read_h5ad(file)
        obsname_dict[file] = list(adata.obs_names)
        del adata # free up memory for large files
    return obsname_dict


def check_label_overlap(obsname_dict: dict[str, list[str]]):
    """
    Throws an error, if any one list contains duplicates

    Prints:
    total amount of unique labels across all files [1,2,3,4] [4,5,6] [3,4,5] -> [1,2,3,4,5,6] -> 6
    
    amount of unique labels shared across all files [1,2,3,4] [4,5,6] [3,4,5] -> [4] -> 1
    percentage of shared labels of unique labels -> 1/6 = 0.166

    amount of labels present in > 1 file [1,2,3,4] [4,5,6] [3,4,5] -> [3,4,5] -> 3
    percentage of labels present in > 1 file of unique labels -> 3/6 = 0.5
    """



    # check if there are duplicates within each list, if so throw error
    def check_duplicates(obsname_dict: dict[str, list[str]]):
        for filepath, obs_names in obsname_dict.items():
            if len(obs_names) != len(set(obs_names)):
                raise ValueError(f"Duplicate labels found in obs_names of {filepath}")
        print("No duplicates found in within each file.")

    # return amount of all labels across all files
    def get_total_labels(obsname_dict: dict[str, list[str]]):
        total_labels = set()
        for obs_names in obsname_dict.values():
            total_labels.update(set(obs_names))
        total_labels_count = len(total_labels)
        return total_labels_count

    # list of labels shared across all files
    def get_shared_labels(obsname_dict: dict[str, list[str]]):
        it = iter(obsname_dict.values())       
        shared_labels = set(next(it)) # get first obs_names list      
        for obs_names in it:                    
            shared_labels &= set(obs_names) # intersect with rest of lists
        return shared_labels

    # list of labels present in > 1 file
    def get_duplicates(obsname_dict: dict[str, list[str]]):
        seen = Counter()
        for lst in obsname_dict.values():
            seen.update(set(lst))   # convert each list to a set first
        duplicates = {item for item, count in seen.items() if count > 1}
        return duplicates
    

    check_duplicates(obsname_dict)
    total_labels_count = get_total_labels(obsname_dict)
    shared_labels_count = len(get_shared_labels(obsname_dict))
    duplicates = get_duplicates(obsname_dict)


    print(f"Total labels: {total_labels_count}")
    print(f"Labels shared across all files: {shared_labels_count}")
    print(f"Percentage of shared labels: {shared_labels_count / total_labels_count:.2f}")
    print(f"Labels present in > 1 file: {len(duplicates)}")
    print(f"Percentage of labels present in > 1 file: {len(duplicates) / total_labels_count:.3f}")


def get_differently_named_cells(h5ad_files: list[str]):
    """
    Prints the amount of cells that have different ids across multiple adatas.
    A cell is uniquely identified by its expression (-> MAKE SURE expression is unmodified (normalization, removing genes, etc.))

    e.g.
    adata1: cellA : [1,2,3,4,5,6]
    adata2: cell 1: [1,2,3,4,5,6]

    -> the cell [1,2,3,4,5,6] occurs under 2 names
    """
    
    expr_to_names = defaultdict(set)

    for adata_file in h5ad_files:   # adata_files = list of paths to MTX/h5ad files
        adata = sc.read_h5ad(adata_file)  # load one AnnData at a time
        print(f"loaded {adata_file}")
        print(f"number of cells: {adata.n_obs}")
        
        for cell_id, vec in zip(adata.obs_names, adata.X): # for loop over 2d array loops over 1 row at a time -> zip gives tuple of cell id and its expression vector
            # Convert sparse row to dense
            if hasattr(vec, "toarray"):
                vec = vec.toarray().ravel() # if sparse, turn to array, ravel to make sure it will be 1D
            key = tuple(vec)  # make it hashable by turning 1D array into tuple (dict can only take hashable types as keys)
            expr_to_names[key].add(cell_id)
        
        del adata  # free memory before loading the next file

    # detect expressions with multiple names
    inconsistent = {k: names for k, names in expr_to_names.items() if len(names) > 1}

    print(f"Amount of cells that exist under different names in different adatas: {len(inconsistent)}")
    print(f"Lenght of expr_to_names: {len(expr_to_names)}")

    import json
    export_dict = {str(key): value for key, value in expr_to_names.items()}

    with open("differently_named_cells.json", "w") as f:
        json.dump(export_dict, f, indent=4)


def get_cells_with_diff_expr(file_list: list[str], layer):
    """
    check if any cells, identified by their cell id, occur under different expression vectors in different adatas (RUNs)
    eg cell a might be (1,2,3,4) in run 1 and (4,3,2,1) in run 2, which should not happen and indicates issues with the pipeline (a cell should always be identical to itself, regardless of processing)
    """


    name_to_vec = defaultdict(set)
    inconsistent = {}

    for file in file_list:
        adata = sc.read_h5ad(file)
        extracted_adata = hf.matrix_to_anndata(adata, layer).copy()
        del adata

        # use csr matrix to save memory (matrix_to_anndata likely outputs np.ndarray for adata.X)
        if type (extracted_adata.X) != sparse.csr_matrix:
            print(f"Type of adata.X: {type(extracted_adata.X)}, converting to csr matrix")
            extracted_adata.X = sparse.csr_matrix(extracted_adata.X)

        for cell_id, vec in zip(extracted_adata.obs_names, extracted_adata.X):

            # turn cell id to string so its hashable
            cell_id = str(cell_id) 

            # Convert sparse row to dense and ravel (to make it 1D)
            if hasattr(vec, "toarray") and sparse.issparse(vec):
                vec = vec.toarray().ravel()
            elif type(vec) == np.ndarray:
                vec = vec.ravel()
            assert type(vec) == np.ndarray

            # convert vec to bytes to save memory, then hash to save even more memory
            vec_bytes = vec.tobytes()
            vec_hash = hashlib.sha256(vec_bytes).digest()

            if cell_id not in inconsistent.keys(): # if we already know the cell has diff expr vecs, no need to add it
                name_to_vec[cell_id].add(vec_hash)

            # check if the current cell_id has > 1  unique expression vector
            if len(name_to_vec[cell_id]) > 1:
                inconsistent[cell_id] = True # inconsistent will only have keys for inconsistent cells, all with value True
                del name_to_vec[cell_id] # delete key from dict to save memory

        del extracted_adata

    print(f"Amount of cells that have different expression in different adatas: {len(inconsistent)}")


if __name__ == "__main__":

    """dir_list = [
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/aggregated",
    ]
    file_list = [os.path.join(dir, file) for dir in dir_list for file in os.listdir(dir)]"""

    file_list = [
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/RUN3.5/reduced/reduced_PDAC_ductal_cell.h5ad",
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/RUN3.6/reduced/reduced_PDAC_ductal_cell.h5ad",
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/RUN4/reduced/reduced_PDAC_ductal_cell.h5ad",
    ]

    get_cells_with_diff_expr(file_list, "X_scANVI_corrected")


    # RESULT: (from 3 aggregated files with slightly different preprocessing paramters, which I thought changes cell order -> changes assigned names in concatenation)
    """
    Getting obs names
    Checking label overlap...
    No duplicates found in within each file.
    Total labels: 48997
    Labels shared across all files: 34074
    Percentage of shared labels: 0.70
    Labels present in > 1 file: 42478
    Percentage of labels present in > 1 file: 0.867
    getting differently named cells
    Amount of cells that exist under different names in different adatas: 0 -> naming is consistent, even with oirignal approach, so we can compare cells by obsname
    """