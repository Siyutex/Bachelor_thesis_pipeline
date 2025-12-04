# this script compares which cells (by labels) are assigned to the transition
# state across runs (issue: labelling might be different)

import scanpy as sc
import helper_functions as hf
import os
from collections import Counter


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
    
    from collections import defaultdict
    expr_to_names = defaultdict(set)

    for adata_file in h5ad_files:   # adata_files = list of paths to MTX/h5ad files
        adata = sc.read_h5ad(adata_file)  # load one AnnData at a time
        
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


if __name__ == "__main__":

    file_list = [
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/isolated/isoltated_test0_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected_cancer_state_inferred_tree_is_['transitional'].h5ad",
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/isolated/isoltated_test1_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected_cancer_state_inferred_tree_is_['transitional'].h5ad",
        r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/RUN4/isolated/isolated_PDAC_ductal_cell_HVG_X_is_X_scANVI_corrected_cancer_state_inferred_tree_is_['transitional'].h5ad",
    ]

    print("Getting obs names")
    obsname_dict = load_obsname_dict(file_list)
    print("Checking label overlap...")
    check_label_overlap(obsname_dict)


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