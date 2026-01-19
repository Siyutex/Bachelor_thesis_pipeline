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
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import json
from matplotlib_venn import venn3, venn2
from typing import Literal
from scipy.optimize import curve_fit
from sklearn.metrics import adjusted_rand_score
from itertools import combinations


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


def load_varname_dict(file_list: list[str]):
    """
    takes a list of filepaths to h5ad files
    return a dictionary, key = filepath, values list(adata.var_names)
    var_names should be ENSG IDs to guarantee unique naming for each gene / feature
    """
    # returns a dict, key = origin file, value = obs_names
    varname_dict = {}
    for file in file_list:
        adata = sc.read_h5ad(file)
        varname_dict[file] = list(adata.var_names)
        del adata # free up memory for large files
    return varname_dict


def check_label_overlap(label_dict: dict[str, list[str]], verbose: bool = False):
    """
    Throws an error, if any one list (value in label_dict) contains duplicates.
    Can be used for obs_names or var_names.

    Prints:
    total amount of unique labels across all files [1,2,3,4] [4,5,6] [3,4,5] -> [1,2,3,4,5,6] -> 6
    
    amount of unique labels shared across all files [1,2,3,4] [4,5,6] [3,4,5] -> [4] -> 1
    percentage of shared labels of unique labels -> 1/6 = 0.166

    amount of labels present in > 1 file [1,2,3,4] [4,5,6] [3,4,5] -> [3,4,5] -> 3
    percentage of labels present in > 1 file of unique labels -> 3/6 = 0.5

    Returns:
    percentage of labels (obsnames / varnames) shared across all files
    """
    vprint = hf.make_vprint(verbose)


    # check if there are duplicates within each list, if so throw error
    def check_duplicates(label_dict: dict[str, list[str]]):
        for filepath, labels in label_dict.items():
            if len(labels) != len(set(labels)):
                raise ValueError(f"Duplicate labels found in file {filepath}")
        vprint("No duplicates found in within each file.")

    # return amount of all labels across all files
    def get_total_labels(label_dict: dict[str, list[str]]):
        total_labels = set()
        for labels in label_dict.values():
            total_labels.update(set(labels))
        total_labels_count = len(total_labels)
        return total_labels_count

    # list of labels shared across all files
    def get_shared_labels(label_dict: dict[str, list[str]]):
        it = iter(label_dict.values())       
        shared_labels = set(next(it)) # get first obs_names list      
        for labels in it:                    
            shared_labels &= set(labels) # intersect with rest of lists
        return shared_labels

    # list of labels present in > 1 file
    def get_duplicates(label_dict: dict[str, list[str]]):
        seen = Counter()
        for lst in label_dict.values():
            seen.update(set(lst))   # convert each list to a set first
        duplicates = {item for item, count in seen.items() if count > 1}
        return duplicates
    

    check_duplicates(label_dict)
    total_labels_count = get_total_labels(label_dict)
    shared_labels_count = len(get_shared_labels(label_dict))
    duplicates = get_duplicates(label_dict)


    vprint(f"Total labels: {total_labels_count}")
    vprint(f"Labels shared across all files: {shared_labels_count}")
    vprint(f"Percentage of shared labels: {shared_labels_count / total_labels_count:.2f}")
    vprint(f"Labels present in > 1 file: {len(duplicates)}")
    vprint(f"Percentage of labels present in > 1 file: {len(duplicates) / total_labels_count:.3f}")

    return shared_labels_count / total_labels_count


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
    
    returns dict, key = cell id, value = true; only contains cells with inconsistent expression across runs
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
    return inconsistent


def get_cells_with_diff_obs(file_list: list[str]) -> dict[True]:
    """
    Check whether any cells, identified by their cell ID, have inconsistent
    tuples of obs annotations across different runs.

    For each adata in file_list, this function extracts for every cell the tuple
    of all obs annotations (i.e. the row in adata.obs). If the same cell has
    different annotation tuples across runs, this indicates nondeterminism or
    issues in the pipeline.

    Returns
    -------
    dict
        Keys are cell IDs (as strings). Values are True.
        Only contains cells whose obs annotation tuples differ across runs.
    """

    name_to_tuple = defaultdict(set)
    inconsistent = {}

    for file in file_list:
        adata = sc.read_h5ad(file)

        # Extract a stable column order once per file
        obs_cols = list(adata.obs.columns)
        obs_df = adata.obs[obs_cols]

        for cell_id, row in zip(adata.obs_names, obs_df.itertuples(index=False, name=None)):
            # Convert cell_id to string for hashing consistency
            cell_id = str(cell_id)

            # row is already a tuple, but its contents may be non-hashable (e.g., arrays)
            # Convert entries to stable byte representation via repr
            normalized = tuple(repr(x) for x in row)

            # Hash tuple to reduce memory footprint
            tuple_bytes = repr(normalized).encode("utf-8")
            tuple_hash = hashlib.sha256(tuple_bytes).digest()

            if cell_id not in inconsistent:
                name_to_tuple[cell_id].add(tuple_hash)

            # Check uniqueness
            if len(name_to_tuple[cell_id]) > 1:
                inconsistent[cell_id] = True
                del name_to_tuple[cell_id]  # free memory

        del adata

    return inconsistent


def get_genes_with_diff_var(file_list: list[str]) -> dict[True]:
    """
    Check whether any genes, identified by their ENSG ID (adata.var_names), have inconsistent
    tuples of var annotations across different runs.

    For each adata in file_list, this function extracts for every gene the tuple
    of all var annotations (i.e. the row in adata.var). If the same gene has
    different annotation tuples across runs, this indicates nondeterminism or
    issues in the pipeline.

    Returns
    -------
    dict
        Keys are ENSG IDs (as strings). Values are True.
        Only contains cells whose var annotation tuples differ across runs.
    """

    name_to_tuple = defaultdict(set)
    inconsistent = {}

    for file in file_list:
        adata = sc.read_h5ad(file)

        # Extract a stable column order once per file
        var_cols = list(adata.var.columns)
        var_df = adata.var[var_cols]

        for gene_id, row in zip(adata.var_names, var_df.itertuples(index=False, name=None)):
            # Convert cell_id to string for hashing consistency
            gene_id = str(gene_id)

            # row is already a tuple, but its contents may be non-hashable (e.g., arrays)
            # Convert entries to stable byte representation via repr
            normalized = tuple(repr(x) for x in row)

            # Hash tuple to reduce memory footprint
            tuple_bytes = repr(normalized).encode("utf-8")
            tuple_hash = hashlib.sha256(tuple_bytes).digest()

            if gene_id not in inconsistent:
                name_to_tuple[gene_id].add(tuple_hash)

            # Check uniqueness
            if len(name_to_tuple[gene_id]) > 1:
                inconsistent[gene_id] = True
                del name_to_tuple[gene_id]  # free memory

        del adata

    return inconsistent


def assert_obsnames_match(dir: str, n_outputs: int = 1, file_name: str = None):
    """
    Assert that all files in a directory sharing the same base `file_name` contain
    identical cell IDs.

    This is useful when a pipeline step produces multiple outputs and you want to
    ensure that corresponding files across runs align (e.g., `run0_*` vs `run1_*`).

    Examples
    --------
    If `file_name="PDAC_cancerous"` and the directory contains:

        run0_PDAC_cancerous_1.h5ad
        run1_PDAC_cancerous_1.h5ad
        run0_PDAC_cancerous_2.h5ad
        run1_PDAC_cancerous_2.h5ad

    then the function checks that:
    - `run0_PDAC_cancerous_1` has the same cell IDs as `run1_PDAC_cancerous_1`
    - `run0_PDAC_cancerous_2` has the same cell IDs as `run1_PDAC_cancerous_2`
    and so on.

    If `file_name` is omitted, the function asserts that *all* files in the
    directory share the same cell IDs. This is useful when the pipeline step
    produces only a single output file per run.

    Parameters
    ----------
    dir : str
        Directory containing the files to compare.
    n_outputs : int, optional
        Number of output files produced per pipeline run. Default is 1.
    file_name : str, optional
        Common name prefix of the files to compare. Assumes file structure of the
        form ``{file_name}_{index}``. For each index, all runs are compared:
        e.g., ``run0_A_0`` vs ``run1_A_0`` vs ``run2_A_0``, etc.
    """


    # create list of all files in directory
    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".h5ad" in file)]


    if file_name != None: # if file_name is given, compare 
        for dataset in range(n_outputs):
            local_file_list = [file for file in file_list if f"{file_name}_{dataset}" in os.path.basename(file)]
            obs_name_dict = load_obsname_dict(local_file_list)
            overlap_percentage = check_label_overlap(obs_name_dict, verbose=False)
            assert overlap_percentage == 1
            print(f"Check passed for {file_name}_{dataset}, cell IDs match across all runs.")
    else:
        obs_name_dict = load_obsname_dict(file_list)
        overlap_percentage = check_label_overlap(obs_name_dict, verbose=False)
        assert overlap_percentage == 1
        print("Check passed, cell IDs match across all files.")


def assert_no_inconsistent_expression(dir: str, layer: str = "X", n_outputs: int = 1, file_name: str = None):
    """
    Check whether the expression vectors of all cells are identical across runs.

    If a pipeline step produces multiple output files per run, you may pass only a
    representative sample, as this function can be slow. For example, if the step
    outputs:

        run0_cell_type_annotated_PDAC_cancerous_0.h5ad
        run0_cell_type_annotated_PDAC_cancerous_1.h5ad
        ...
        run1_cell_type_annotated_PDAC_cancerous_0.h5ad
        run1_cell_type_annotated_PDAC_cancerous_1.h5ad
        ...
        run2_cell_type_annotated_PDAC_cancerous_0.h5ad
        run2_cell_type_annotated_PDAC_cancerous_1.h5ad

    then you only need to pass:

        run0_cell_type_annotated_PDAC_cancerous_0.h5ad
        run1_cell_type_annotated_PDAC_cancerous_0.h5ad
        run2_cell_type_annotated_PDAC_cancerous_0.h5ad

    If a single set of corresponding files matches across runs, this should be
    sufficient to verify determinism.

    Parameters
    ----------
    dir : str
        Directory containing the files to compare.
    n_outputs : int, optional
        Number of output files to compare per run. Default is 1.
        Since this function is slow, using ``n_outputs > 1`` is not recommended.
    file_name : str, optional
        Common name prefix of the files to compare. Assumes a file structure of the
        form ``{file_name}_{index}``. For each index (0 to ``n_outputs - 1``),
        all runs are compared—for example: ``run0_A_0`` vs. ``run1_A_0`` vs.
        ``run2_A_0``.
    layer: str, optional
        Layer to compare. Default is "X"
        Cant take "X", any adata.layers matrix, or any adata.obsm matrix
    """

    # create list of all files in directory
    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".h5ad" in file)]

    if file_name != None: # if file_name is given, compare 
        for dataset in range(n_outputs):
            local_file_list = [file for file in file_list if f"{file_name}_{dataset}" in os.path.basename(file)]
            inconsistent_dict = get_cells_with_diff_expr(local_file_list, layer)
            assert len(inconsistent_dict) == 0
            print(f"Check passed for {file_name}_{dataset}, expression matches across all runs for all cell IDs in layer {layer}.")
    else:
        inconsistent_dict = get_cells_with_diff_expr(file_list, layer)
        assert len(inconsistent_dict) == 0
        print(f"Check passed, expression matches across all runs for all cell IDs in layer {layer}.")


def assert_no_inconsistent_obs_annotations(dir: str, n_outputs: int = 1, file_name: str = None):
    """
    Check whether the obs annotations of all cells are identical across runs.

    If a pipeline step produces multiple output files per run, you may pass only a
    representative sample, as this function can be slow. For example, if the step
    outputs:

        run0_cell_type_annotated_PDAC_cancerous_0.h5ad
        run0_cell_type_annotated_PDAC_cancerous_1.h5ad
        ...
        run1_cell_type_annotated_PDAC_cancerous_0.h5ad
        run1_cell_type_annotated_PDAC_cancerous_1.h5ad
        ...
        run2_cell_type_annotated_PDAC_cancerous_0.h5ad
        run2_cell_type_annotated_PDAC_cancerous_1.h5ad

    then you only need to pass:

        run0_cell_type_annotated_PDAC_cancerous_0.h5ad
        run1_cell_type_annotated_PDAC_cancerous_0.h5ad
        run2_cell_type_annotated_PDAC_cancerous_0.h5ad

    If a single set of corresponding files matches across runs, this should be
    sufficient to verify determinism.

    Parameters
    ----------
    dir : str
        Directory containing the files to compare.
    n_outputs : int, optional
        Number of output files to compare per run. Default is 1.
        Since this function is slow, using ``n_outputs > 1`` is not recommended.
    file_name : str, optional
        Common name prefix of the files to compare. Assumes a file structure of the
        form ``{file_name}_{index}``. For each index (0 to ``n_outputs - 1``),
        all runs are compared—for example: ``run0_A_0`` vs. ``run1_A_0`` vs.
        ``run2_A_0``.
    """

    # create list of all files in directory
    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".h5ad" in file)]

    if file_name != None: # if file_name is given, compare 
        for dataset in range(n_outputs):
            local_file_list = [file for file in file_list if f"{file_name}_{dataset}" in os.path.basename(file)]
            inconsistent_dict = get_cells_with_diff_obs(local_file_list)
            assert len(inconsistent_dict) == 0
            print(f"Check passed for {file_name}_{dataset}, obs annotations match across all runs.")
    else:
        inconsistent_dict = get_cells_with_diff_obs(file_list)
        assert len(inconsistent_dict) == 0
        print("Check passed, obs annotations match across all files.")


def assert_varnames_match(dir: str, n_outputs: int = 1, file_name: str = None):
    """
    Assert that all files in a directory sharing the same base `file_name` contain
    identical adata.var_names.

    This is useful when a pipeline step produces multiple outputs and you want to
    ensure that corresponding files across runs align (e.g., `run0_*` vs `run1_*`).

    Examples
    --------
    If `file_name="PDAC_cancerous"` and the directory contains:

        run0_PDAC_cancerous_1.h5ad
        run1_PDAC_cancerous_1.h5ad
        run0_PDAC_cancerous_2.h5ad
        run1_PDAC_cancerous_2.h5ad

    then the function checks that:
    - `run0_PDAC_cancerous_1` has the same cell IDs as `run1_PDAC_cancerous_1`
    - `run0_PDAC_cancerous_2` has the same cell IDs as `run1_PDAC_cancerous_2`
    and so on.

    If `file_name` is omitted, the function asserts that *all* files in the
    directory share the same cell IDs. This is useful when the pipeline step
    produces only a single output file per run.

    Parameters
    ----------
    dir : str
        Directory containing the files to compare.
    n_outputs : int, optional
        Number of output files produced per pipeline run. Default is 1.
    file_name : str, optional
        Common name prefix of the files to compare. Assumes file structure of the
        form ``{file_name}_{index}``. For each index, all runs are compared:
        e.g., ``run0_A_0`` vs ``run1_A_0`` vs ``run2_A_0``, etc.
    """


    # create list of all files in directory
    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".h5ad" in file)]


    if file_name != None: # if file_name is given, compare 
        for dataset in range(n_outputs):
            local_file_list = [file for file in file_list if f"{file_name}_{dataset}" in os.path.basename(file)]
            var_name_dict = load_varname_dict(local_file_list)
            overlap_percentage = check_label_overlap(var_name_dict, verbose=False)
            assert overlap_percentage == 1
            print(f"Check passed for {file_name}_{dataset}, var names match across all runs.")
    else:
        var_name_dict = load_varname_dict(file_list)
        overlap_percentage = check_label_overlap(var_name_dict, verbose=False)
        assert overlap_percentage == 1
        print("Check passed, var names match across all files.")


def assert_no_inconsistent_var_annotations(dir: str, n_outputs: int = 1, file_name: str = None):
    """
    Check whether the var annotations of all genes are identical across runs.

    If a pipeline step produces multiple output files per run, you may pass only a
    representative sample, as this function can be slow. For example, if the step
    outputs:

        run0_cell_type_annotated_PDAC_cancerous_0.h5ad
        run0_cell_type_annotated_PDAC_cancerous_1.h5ad
        ...
        run1_cell_type_annotated_PDAC_cancerous_0.h5ad
        run1_cell_type_annotated_PDAC_cancerous_1.h5ad
        ...
        run2_cell_type_annotated_PDAC_cancerous_0.h5ad
        run2_cell_type_annotated_PDAC_cancerous_1.h5ad

    then you only need to pass:

        run0_cell_type_annotated_PDAC_cancerous_0.h5ad
        run1_cell_type_annotated_PDAC_cancerous_0.h5ad
        run2_cell_type_annotated_PDAC_cancerous_0.h5ad

    If a single set of corresponding files matches across runs, this should be
    sufficient to verify determinism.

    Parameters
    ----------
    dir : str
        Directory containing the files to compare.
    n_outputs : int, optional
        Number of output files to compare per run. Default is 1.
        Since this function is slow, using ``n_outputs > 1`` is not recommended.
    file_name : str, optional
        Common name prefix of the files to compare. Assumes a file structure of the
        form ``{file_name}_{index}``. For each index (0 to ``n_outputs - 1``),
        all runs are compared—for example: ``run0_A_0`` vs. ``run1_A_0`` vs.
        ``run2_A_0``.
    """

    # create list of all files in directory
    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".h5ad" in file)]

    if file_name != None: # if file_name is given, compare 
        for dataset in range(n_outputs):
            local_file_list = [file for file in file_list if f"{file_name}_{dataset}" in os.path.basename(file)]
            inconsistent_dict = get_genes_with_diff_var(local_file_list)
            assert len(inconsistent_dict) == 0
            print(f"Check passed for {file_name}_{dataset}, var annotations match across all runs.")
    else:
        inconsistent_dict = get_genes_with_diff_var(file_list)
        assert len(inconsistent_dict) == 0
        print("Check passed, var annotations match across all files.")


def run_all_h5ad_checks(dir: str, n_outputs: int = 1, file_name: str = None, layer="X"):
    """
    Run alle assertion functions:
    - assert_obsnames_match (uses passed n_outputs)
    - assert_no_inconsistent_expression (only compares a representative sample n_outputs = 1)
    - assert_no_inconsistent_obs_annotations (only compares a representative sample n_outputs = 1)
    - assert_varnames_match (uses passed n_outputs)
    - assert_no_inconsistent_var_annotations (only compares a representative sample n_outputs = 1)

    If any of them fail, print a message indicating which check failed.

    Parameters
    ----------
    dir : str
        Directory containing the files to compare.
    n_outputs : int, optional
        Number of output files to compare per run. Default is 1.
    file_name : str, optional
        Common name prefix of the files to compare. Assumes a file structure of the
        form ``{file_name}_{index}``. For each index (0 to ``n_outputs - 1``),
        all runs are compared—for example: ``run0_A_0`` vs. ``run1_A_0`` vs.
        ``run2_A_0``.
        Set to None to compare all files in the directory (overrides n_outputs).
    layer : str, optional
        Layer to compare expression of, for assert_no_inconsistent_expression. Default is "X".
    """

    try:
        assert_obsnames_match(dir, n_outputs=n_outputs, file_name=file_name)
    except AssertionError:
        print("\n CHECK FAILED, obs names do not match across runs. \n")
        pass
    try:
        assert_no_inconsistent_expression(dir, layer=layer, n_outputs=1, file_name=file_name)
        pass
    except AssertionError:
        print("\n CHECK FAILED, expression matrices do not match across runs. \n")
        pass
    try:
        assert_no_inconsistent_obs_annotations(dir, n_outputs=1, file_name=file_name)
        pass
    except AssertionError:
        print("\n CHECK FAILED, obs annotations do not match across runs. \n")
        pass
    try:
        assert_varnames_match(dir, n_outputs=n_outputs, file_name=file_name)
    except AssertionError:
        print("\n CHECK FAILED, var names do not match across runs. \n")
        pass
    try:
        assert_no_inconsistent_var_annotations(dir, n_outputs=1, file_name=file_name)
    except AssertionError:
        print("\n CHECK FAILED, var annotations do not match across runs. \n")
        pass


def assert_tree_equivalence(dir: str):
    """
    Assert that the trees in nwk files in the directory directory are identical.
    Useful for comparing non h5ad files, like nwk.
    """
    import skbio

    def strip_distances(tree):
        for node in tree.traverse():
            node.length = None
        return tree

    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".nwk" in file)]
    
    # get tree strings without distances in a list
    str_list = []
    for file in file_list:
        tree = skbio.TreeNode.read(file)
        tree = strip_distances(tree)
        str_list.append(str(tree))

    # get list of hashes of trees
    hash_list = []
    for tree_string in str_list:
        tree_bits = tree_string.encode("utf-8")
        tree_hash = hashlib.sha256(tree_bits).digest()
        hash_list.append(tree_hash)

    # compare hashes
    inconsistent = False
    for i, hash in enumerate(hash_list):
        print(f"hash {i} == hash 0: {hash == hash_list[0]}")
        if hash != hash_list[0]:
            inconsistent = True
            inconsistent_tree_index = i

    # evauluate results
    if inconsistent == False: # if hashes are the same, print message
        print("Check passed, all trees are equivalent.")
        return
    else: # else find differences
        for i in range(max(len(str_list[0]), len(str_list[inconsistent_tree_index]))):
            if str_list[0][i] != str_list[inconsistent_tree_index][i]:
                char = str_list[0][i]
                char_zero = str_list[inconsistent_tree_index][i]
                print(f"{char_zero} != {char}")
                print(f"Surrounding text: {str_list[inconsistent_tree_index][i-10:i+10]}")
                return


def evaluate_TS_consistency(dir: str):
    """
    Take a number of subsampled files with transition states annotated.
    Check consistency of TSs across runs:
        - amount (abs / rel) of cells always labelled as TS
        - amount (abs / rel) of cells always not labelled as TS
        - amount (abs / rel) of cells that are labelled as TS at at least one run, but not in all runs
    """  

    # create dataframe that holds cell states for all cells for each run (or "not_sampled" if the cell is not present in the run)
    def make_df_plots(state_df: pd.DataFrame):
        """
        take dataframe with rows = cells, columns = runs, values = state

        make a heatmap of probability of each state for each cell (probability calculated from relative frequency of each state across runs where the cell was sampled)
        make a colormap of cell states across runs for all cells
        """
        
        def plot_probability_df(df):
            """
            input is a dataframe with rows = cells, columns = runs, values = state / not_sampled

            make and save a heatmap of P(state|sampled) for states: cancerous, non_cancerous, transitional
            if states are consistent across runs, each cell should have 1 for one state and 0 for all other states
            """
            
            def get_state_probabilities(df, state: str):
                """
                returns a vector of P(state|sampled) with a value for each cell across all runs
                """
                state_counts = df.eq(state).sum(axis=1)
                not_sampled_counts = df.eq("not_sampled").sum(axis=1)
                state_probabilities = state_counts / (df.shape[1] - not_sampled_counts)

                return state_probabilities

            df = pd.DataFrame({
                "cancerous": get_state_probabilities(df, "cancerous"),
                "non_cancerous": get_state_probabilities(df, "non_cancerous"),
                "transitional": get_state_probabilities(df, "transitional"),
            })

            sns.heatmap(df, annot=False, cmap="coolwarm", vmax=1, vmin=0)
            plt.savefig(os.path.join(dir, "State_probabilities.png"))

        def plot_state_df(df):
            """
            input is a dataframe with rows = cells, columns = runs, values = state / not_sampled

            make and save colormap of states across runs for each cell
            eg cell 1: {run1: cancerous -> blue, run2: non_canerous -> orange, run3: transitional -> green}
            """

            state_colors = { # hex colors
                "cancerous": "#1472bf",
                "non_cancerous": "#ff7f00",
                "transitional": "#42ab3e",
                "unassigned": "#000000",
                "not_sampled": "#FFFFFF"
            }

            # mapping from state to int, sns needs numeric values
            state_to_int = {state: i for i, state in enumerate(state_colors.keys())}
            heatmap_data = df.replace(state_to_int)

            # color map (list of hex colors in same order as state_to_int sorted by int value)
            colors = [state_colors[state] for state, i in sorted(state_to_int.items(), key=lambda x: x[1])]
            cmap = ListedColormap(colors)

            sns.heatmap(
                heatmap_data,
                annot=False, # do not show orignal labels on data, we have too many cells for that
                cmap=cmap,
                cbar=False,  # do not show color bar
            )
            plt.ylabel("Cell ID")
            plt.xlabel("Run")
            plt.title("Cell states across runs")
            plt.savefig(os.path.join(dir, "State_colormap.png"))
                    

        print("plot probability_df")
        plot_probability_df(state_df)
        print("plot state_df")
        plot_state_df(state_df)


    # get list of h5ad files corresponding to runs
    file_list = [os.path.join(dir, file) for file in os.listdir(dir) if (os.path.isfile(os.path.join(dir, file)) and ".h5ad" in file)]

    # get dict of run: [(cell ID, state), ... ]
    print("cerating run_dict")
    run_dict = {}
    for file in file_list:
        ad = sc.read_h5ad(file)
        run_dict[file] = list(zip(ad.obs_names, ad.obs["cancer_state_inferred_tree"]))
    # get dict for run: [cell ID, ...] if that cell ID corresponds to a TS cell
    print("creating TS_cells_per_run")
    TS_cells_per_run = {
        run: [cell_id for (cell_id, state) in cells if state == "transitional"]
        for run, cells in run_dict.items()
    }
    # create dataframe with columns = runs, rows = cells, values = state
    print("creating state_df")
    dfs = []
    for run_id, cells in run_dict.items():
        df = pd.DataFrame(cells, columns=["cell_id", run_id]) # use cell ID as column for now, needed for df merging, run_id = in this run this cell had that value
        dfs.append(df)
    state_df = dfs[0]
    for df in dfs[1:]:
        state_df = pd.merge(state_df, df, on="cell_id", how="outer") # merge dfs so each cell_Id (row) always gets the data corresponding to it
    state_df = state_df.set_index("cell_id") # set cell ID as index and remove the cell ID column
    state_df = state_df.fillna("not_sampled") # if a cell was not present in a run then label it as "not_sampled"

    # make heatmap of state probability per cell and colormap of state per cell across runs
    make_df_plots(state_df) 


    # get list of all cell IDs that occur across runs
    all_cell_IDs = set()
    for cells in run_dict.values():
        all_cell_IDs.update([cells[i][0] for i in range(len(cells))])

    # get list of all TS cells that occur across all runs in which it was sampled
    ts_counts = state_df.eq("transitional").sum(axis=1) # count number of times each cell was labelled as TS
    not_sampled_counts = state_df.eq("not_sampled").sum(axis=1)
    consistent_ts_cells_mask = ts_counts == state_df.shape[1] - not_sampled_counts # boolean mask
    consistent_ts_cells = set(state_df[consistent_ts_cells_mask].index)

    # get list of all TS cells that occur at least once, but not always (inconsistently labelled cells)
    inconsistent_TS_cells = set.union(*(set(value) for value in TS_cells_per_run.values())) # all TS cells that ever occured (*(generator) unpacks the generator and passes each generated expression as a seperate argument)
    inconsistent_TS_cells -= consistent_ts_cells # remove cells that are consistently labelled as TS

    # get list of all cells that never get labelled as TS
    consistent_non_TS_cells = all_cell_IDs - (consistent_ts_cells | inconsistent_TS_cells) 

    print(f"Number of total unique cells: {len(all_cell_IDs)}")
    print(f"\nNumber of consistently labelled transitional cells:{len(consistent_ts_cells)}")
    print(f"Realtive amount of consistently labelled transitional cells: {len(consistent_ts_cells)/len(all_cell_IDs)}")
    print(f"\nNumber of inconsistently labelled transitional cells: {len(inconsistent_TS_cells)}")
    print(f"Relative amount of inconsistently labelled transitional cells: {len(inconsistent_TS_cells)/len(all_cell_IDs)}")
    print(f"\nNumber of consistently labelled non-transitional cells: {len(consistent_non_TS_cells)}")
    print(f"Relative amount of consistently labelled non-transitional cells: {len(consistent_non_TS_cells)/len(all_cell_IDs)}")

    




def evaluate_set_consistency(directory_path, set_type: Literal["edges", "var_names", "obs_names"] ):
    """
    For set_type = "edges":
    Take 3 json files with edges from differnt GRN edge inference runs.
    Compute pairwise and global jaccard similarity + create venn diagram of edge overlap.

    For set_type = "var_names":
    Take 3 h5ad files with variable names from different pipeline runs.
    Compute pairwise and global jaccard similarity + create venn diagram of variable name overlap.

    For set_type = "obs_names":
    Take 3 h5ad files with observation names from different pipeline runs.
    Compute pairwise and global jaccard similarity + create venn diagram of observation name overlap.
    """


    def load_edges_from_json(filepath):
        """
        Parses the JSON and flattens it into a set of directed edge tuples.
        Example: {"Source": ["T1", "T2"]} -> {("Source", "T1"), ("Source", "T2")}
        """
        edges = set()
        with open(filepath, 'r') as f:
            data = json.load(f)
            for source, targets in data.items():
                for target in targets:
                    edges.add((source, target))
        return edges
    
    def load_var_names_from_h5ad(filepath):
        ad = sc.read_h5ad(filepath)
        return set(ad.var_names)
    
    def load_obs_names_from_h5ad(filepath):
        ad = sc.read_h5ad(filepath)
        return set(ad.obs_names)

    def calculate_metrics(sets_list, filenames):
        # Pairwise Jaccard for reference
        print("--- Pairwise Jaccard Similarity ---")
        for i in range(len(sets_list)):
            for j in range(i + 1, len(sets_list)):
                s1, s2 = sets_list[i], sets_list[j]
                intersection = len(s1.intersection(s2))
                union = len(s1.union(s2))
                jaccard = intersection / union if union > 0 else 0
                print(f"Run_{i} vs Run_{j}: {jaccard:.4f}")

        # Global Jaccard: (A & B & C) / (A | B | C)
        global_intersection = set.intersection(*sets_list)
        global_union = set.union(*sets_list)
        global_jaccard = len(global_intersection) / len(global_union) if global_union else 0
        
        print("\n--- Global Metrics ---")
        print(f"Global Jaccard (Intersection of all / Union of all): {global_jaccard:.4f}")
        return global_jaccard

    def plot_grn_venn(sets_list, filenames):
        plt.figure(figsize=(10, 8))
        # Create the Venn diagram
        if len(sets_list) == 3:
            v = venn3(sets_list, set_labels=('Run 1', 'Run 2', 'Run 3'))
        elif len(sets_list) == 2:
            v = venn2(sets_list, set_labels=('Run 1', 'Run 2'))
        else:
            raise ValueError("Unsupported number of sets for Venn diagram.")
        
        
        plt.savefig(os.path.join(directory_path, "venn.png"))


    # Get first 3 json or h5ad files
    files = [f for f in os.listdir(directory_path) if f.endswith(".json") or f.endswith(".h5ad")][:3]
    
    if len(files) == 2:
        print(f"Found {len(files)} files. Creating 2 - set Venn diagram.")
    elif len(files) == 3:
        print(f"Found {len(files)} files. Creating 3 - set Venn diagram.")
    else:
        print(f"Found {len(files)} files. Cannot create Venn diagram.")
        return # leave function

    sets = []
    for file in files:
        full_path = os.path.join(directory_path, file)
        if set_type == "var_names":
            sets.append(load_var_names_from_h5ad(full_path))
        elif set_type == "obs_names":
            sets.append(load_obs_names_from_h5ad(full_path))
        elif set_type == "edges":
            sets.append(load_edges_from_json(full_path))
        print(f"Loaded {len(sets[-1])} elements from {file}")
    print("\n")

    calculate_metrics(sets, files)
    plot_grn_venn(sets, files)


def find_consistency_limit(dir: str, set_type: Literal["edges", "var_names", "obs_names"] ):


    def load_edges_from_json(filepath):
        """
        Parses the JSON and flattens it into a set of directed edge tuples.
        Example: {"Source": ["T1", "T2"]} -> {("Source", "T1"), ("Source", "T2")}
        """
        edges = set()
        with open(filepath, 'r') as f:
            data = json.load(f)
            for source, targets in data.items():
                for target in targets:
                    edges.add((source, target))
        return edges
    
    def load_var_names_from_h5ad(filepath):
        ad = sc.read_h5ad(filepath)
        return set(ad.var_names)
    
    def load_obs_names_from_h5ad(filepath):
        ad = sc.read_h5ad(filepath)
        return set(ad.obs_names)
    
    files = [f for f in os.listdir(dir) if f.endswith(".json") or f.endswith(".h5ad")]
    sets = []
    for file in files:
        full_path = os.path.join(dir, file)
        if set_type == "var_names":
            sets.append(load_var_names_from_h5ad(full_path))
        elif set_type == "obs_names":
            sets.append(load_obs_names_from_h5ad(full_path))
        elif set_type == "edges":
            sets.append(load_edges_from_json(full_path))
        print(f"Loaded {len(sets[-1])} elements from {file}")

    def avg_pairwise_jaccard(set_list):
        jaccard_list = []
        for i in range(len(set_list)):
            for j in range(i + 1, len(set_list)):
                s1, s2 = set_list[i], set_list[j]
                intersection = len(s1.intersection(s2))
                union = len(s1.union(s2))
                jaccard = intersection / union if union > 0 else 0
                jaccard_list.append(jaccard)
        
        avg = sum(jaccard_list) / len(jaccard_list)
        return avg
    

    def fit_jaccard_limit(avg_jaccards):
        """
        Fits a growth function to average Jaccard scores to find the limit L.
        Model: f(n) = L - a * exp(-b * n)
        """
        n_values = np.arange(2, len(avg_jaccards) + 2)
        y_values = np.array(avg_jaccards)

        # Growth model where L is the upper asymptote
        def growth_model(n, L, a, b):
            return L - a * np.exp(-b * n)

        # Initial guesses:
        # L: slightly higher than the max observed value
        # a: the difference between the limit and the starting point
        # b: a small growth rate
        p0 = [y_values[-1] + 0.05, y_values[-1] - y_values[0], 0.1]
        
        # Constraints: L must be between 0 and 1
        bounds = (0, [1.0, 1.0, np.inf])
        
        try:
            params, _ = curve_fit(growth_model, n_values, y_values, p0=p0, bounds=bounds)
            L, a, b = params
        except Exception as e:
            print(f"Fitting failed: {e}")
            return None

        # --- Visualization ---
        plt.figure(figsize=(10, 6))
        plt.scatter(n_values, y_values, color='red', label='Observed Avg Jaccard')
        
        # Generate curve
        n_smooth = np.linspace(2, len(avg_jaccards) + 5, 100)
        plt.plot(n_smooth, growth_model(n_smooth, L, a, b), 'b--', 
                label=f'Fit Curve (Limit L ≈ {L:.4f})')
        
        plt.axhline(y=L, color='green', linestyle=':', label=f'Asymptote (L={L:.4f})')
        plt.title('Consistency Limit: Average Pairwise Jaccard Convergence')
        plt.xlabel('Number of Sets (n)')
        plt.ylabel('Average Jaccard')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(dir, "jaccard_convergence.png"))

        return L

    # track avg pariwise jaccard from 1 to n sets
    avg_jaccard_list = []
    for i in range(1,len(sets)): # iterate from 1 to last index of list
        avg_jaccard_list.append(avg_pairwise_jaccard(sets[:i+1])) # compare the first i+1 sets (but does not break at end)

    # find limit
    limit = fit_jaccard_limit(avg_jaccard_list)
    print(f"Limit: {limit}")
    

def get_pairwise_ari(directory, resolution=0.5, n_neighbors=15, n_comps=50):
    """
    ARI = adjusted rand index

    ARI = 1 means that a particular cell gets assigned in a cluster with the same other cells across subsampling runs
    ARI = 0 means that cells get assigned to cluster randomly over subsampling runs

    This can be used to evauate large scale topological consistency of a dataset (eg Shin et al data might be less noisy -> subsampling does not affect large scale clustering as much -> more consistent tree -> clades -> TS)

    It basically works by getting a vector for each of 2 runs where the index = the cell and the value = what cluster it is in in that run
    then it compares for all pairs of cells (n_obs over 2): "If cell A and cell B are in the same cluster (name of cluster does not matter) in run 1, or they also in the same cluster in run 2?"
    possible results for this (per pair of cells are): in same cluster in both runs, never in same cluster, once serperated and once together
    in the first 2 cases, we can say that the pair of cells gets consistently clustered together or apart -> agreement

    We then check for this clustering agreement for each pair of cells across the 2 runs
    then RI = number of agreements / number of cell pairs 
    1 = perfect agreement between 2 runs
    0 = random assignment to clusters

    ARI improves this by Adjusting for random noise (randomly distributed cells can also happen to be together / apart in both runs = agreement)
    so ARI tells us how consistent the clustering (=global topology) is compared to random noise

    We then compute this ARI for all pairs of runs (pairwise ARI) (because ARI over all runs would be super super harsh (almost no cells would consistenly group together over 50 runs))
    then we take the average

    isi ai shisgamebre di aisklin

    resoltution = resolution of leiden to use for clustering
    """

    files = [f for f in os.listdir(directory) if f.endswith('.h5ad')]
    
    # 1. Clustering Registry: Store only the final labels to save RAM
    label_registry = {}
    cluster_counts = []
    
    print(f"Processing and clustering {len(files)} files...")
    
    for f in files:
        # Load full file into memory to allow neighbor/leiden computation
        path = os.path.join(directory, f)
        adata = sc.read_h5ad(path)
        adata = hf.matrix_to_anndata(adata, matrix_key="X_scANVI_corrected")
        
        # Recalculate neighbors and Leiden on this specific subsample
        # This ensures the graph is built specifically from the subsampled cells
        sc.pp.pca(adata, n_comps=n_comps, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca")
        sc.tl.leiden(adata, resolution=resolution, key_added='temp_cluster')
        
        # Store the Series (obs_names -> cluster_id)
        label_registry[f] = adata.obs['temp_cluster'].copy()

        # get amount of leiden clusters
        cluster_counts.append(len(adata.obs['temp_cluster'].unique()))
        
        # Clean up to keep memory free for the next file
        del adata
    
    # 2. Fast Pairwise ARI calculation
    results = []
    pairs = list(combinations(files, 2))
    print(f"Calculating ARI for {len(pairs)} pairs...")

    for f1, f2 in pairs:
        labels1 = label_registry[f1]
        labels2 = label_registry[f2]
        
        # Find intersection of cell barcodes
        common_cells = labels1.index.intersection(labels2.index)
        
        # Safety check: ARI requires at least 2 points to form a pair
        if len(common_cells) < 2:
            print(f"Skipping pair {f1} and {f2}: insufficient overlap.")
            continue
            
        # Align labels based on the common cells
        score = adjusted_rand_score(
            labels1.loc[common_cells], 
            labels2.loc[common_cells]
        )
        
        results.append({
            'run_a': f1, 
            'run_b': f2, 
            'ari': score, 
            'n_common': len(common_cells)
        })

    df_results = pd.DataFrame(results)
    
    if not df_results.empty:
        print(f"\nPairwise ARI statistics for {directory}:")
        print(df_results['ari'].describe())
    
    print(f"Cluster statistics"):
    print(f"Mean cluster count: {np.mean(cluster_counts)}")
    print(f"Std cluster count: {np.std(cluster_counts)}")
    print(f"Median cluster count: {np.median(cluster_counts)}")
    
    return df_results



if __name__ == "__main__":

    print("starting script...")

    dir = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/subsampled/PDAC/jaccard_convergency_check"
    dir2 = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/subsampled/Shin"
    get_pairwise_ari(dir)
    get_pairwise_ari(dir2)
    

    