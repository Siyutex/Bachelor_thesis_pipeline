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
    file_list = [os.path.join(dir, file) for file in os.listdir(dir)]


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
    file_list = [os.path.join(dir, file) for file in os.listdir(dir)]

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
    file_list = [os.path.join(dir, file) for file in os.listdir(dir)]

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
    file_list = [os.path.join(dir, file) for file in os.listdir(dir)]


    if file_name != None: # if file_name is given, compare 
        for dataset in range(n_outputs):
            local_file_list = [file for file in file_list if f"{file_name}_{dataset}" in os.path.basename(file)]
            var_name_dict = load_varname_dict(local_file_list)
            overlap_percentage = check_label_overlap(var_name_dict, verbose=False)
            assert overlap_percentage == 1
            print(f"Check passed for {file_name}_{dataset}, varnames match across all runs.")
    else:
        var_name_dict = load_varname_dict(file_list)
        overlap_percentage = check_label_overlap(var_name_dict, verbose=False)
        assert overlap_percentage == 1
        print("Check passed, cell IDs match across all files.")


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
    file_list = [os.path.join(dir, file) for file in os.listdir(dir)]

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


def run_all_checks(dir: str, n_outputs: int = 1, file_name: str = None, layer="X"):
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


if __name__ == "__main__":

    print("starting script...")

    dir = r"/proj/ml_grn/project_julian/Bachelor_thesis_pipeline/Data/output_storage/batch_corrected"
    run_all_checks(dir=dir, n_outputs=1, file_name=None, layer="X_scANVI_corrected")




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