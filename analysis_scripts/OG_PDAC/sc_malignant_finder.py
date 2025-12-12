# uses sc malignant finder to determine malignancy states for each cell 

from scMalignantFinder import classifier
import scanpy as sc
import helper_functions as hf
import warnings
import os


def check_assign_cell_type(adata, cell_type) -> sc.AnnData:

    # limit adata to chosen cell_type if given
    if cell_type is not None and cell_type in adata.obs["cell_type"].values:
        vprint(f"Isolating {cell_type}...")
        internal_adata = adata[adata.obs["cell_type"] == cell_type, :].copy()
    elif cell_type is not None and cell_type not in adata.obs["cell_type"].values:
        raise ValueError(f"Cell type {cell_type} not found in adata.obs['cell_type']")
    elif cell_type is None:
        vprint("Using all cells...") 

    return internal_adata # we need to return the object here because slicing creates a new object, so the adata passed as a paremeter is not modified



def run_scMF(adata, layer) -> sc.AnnData:
    # runs scMF and returns anndata with annotations added

    # define adata and make sure varnames are gene symbols
    internal_adata = hf.matrix_to_anndata(adata, layer)
    if adata.var_names.str.startswith("ENSG").all(): # if var names are ensembl ids
        vprint("Setting varnames to gene symbols (needed for scMF)...")
        internal_adata.var_names = internal_adata.var["gene_symbols"] # unique mapping ensured by aggregation script

    # check for and remvoe duplicate varnames
    vprint("checking for duplicate gene names...")
    if internal_adata.var_names.duplicated().any():
        vprint(f"Duplicate varnames{internal_adata.var_names[internal_adata.var_names.duplicated()]}")
        # if duplicated varnames are < 1% remove them
        if (internal_adata.var_names.duplicated().sum() / len(internal_adata.var_names)) < 0.01:
            vprint("removing duplicate varnames...")
            internal_adata = internal_adata[:, ~internal_adata.var_names.duplicated()].copy()
        else:
            raise ValueError("More than 1% of varnames are duplicates, please validate your input data")

    # check normalization state
    vprint("Checking normalization state...")
    is_normalized = hf.is_normalized(internal_adata)
    if not is_normalized:
        warnings.warn("adata does not seem to be normalized, did you pass the right layer?")

    # set path to pretrained model
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pretrained_model_path = os.path.join(script_dir, "..", "..", "pretrained_scMF_model")
    
    # intialize pretrained model
    model = classifier.scMalignantFinder(
        test_input=internal_adata,
        pretrain_dir=pretrained_model_path, # refernce to dir with pretrained model and feature list
        n_thread=-1, # usually means all cores
        norm_type= not is_normalized # normalize if not already
    )

    # load model
    vprint("loading model")
    model.load()

    # predict
    vprint("predicting malignancy states")
    internal_adata = model.predict().copy() # returns anndata object with malignancy states in obs["scMalignantFinder_prediction"]
    malignancy_state = internal_adata.obs["scMalignantFinder_prediction"]
    malignancy_probabilities = internal_adata.obs["malignancy_probability"] # returns anndata object with malignancy probabilities in obs["scMalignantFinder_probabilities"]
    del internal_adata # free up memory

    # standardize nomenclature
    malignancy_state = malignancy_state.str.replace("Malignant", "cancerous") # pandas series supports replace
    malignancy_state = malignancy_state.str.replace("Normal", "non_cancerous")

    result_adata = adata.copy()
    result_adata.obs["cancer_state_inferred_scMF"] = malignancy_state
    result_adata.obs["malignancy_probability_scMF"] = malignancy_probabilities

    return result_adata


def main(input_data_file, output_data_dir, layer, cell_type):

    # load adata
    print("reading adata")
    adata = sc.read_h5ad(input_data_file)
    vprint(f"Adata summary:\n{adata}")

    # limit cells
    if cell_type:
        print(f"limiting cells to {cell_type}...")
        adata = check_assign_cell_type(adata, cell_type).copy()

    # run scMF
    print("running scMF...")
    adata = run_scMF(adata, layer).copy()

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))





if __name__ == "__main__":

    # import arguments
    input_data_file, output_data_dir, layer, cell_type, verbose = hf.import_cmd_args(5)

    # define vprint
    vprint = hf.make_vprint(verbose)

    main(input_data_file, output_data_dir, layer, cell_type)
