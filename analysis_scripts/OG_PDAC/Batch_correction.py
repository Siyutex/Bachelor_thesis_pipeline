# The input to this script should be a directory with preprocessed h5ad files (only filtering, no normlization, no HVG selection)
# The script aggregates all datasets of a single tissue type (cancer and non cancer) into an anndata
# Then it runs batch correction ber on that anndata using scVI, training a new model if one does not exist for the set parameters
# all of this will happen directly in gene expression space
# The results can also be plotted

# assumes ensebl ids as adata.var_names and gene symbols in adata.var.["gene_symbols"]

import os
import scanpy as sc
import scvi # needed for batch correction
import helper_functions as hf
import anndata
from typing import Tuple


def correct_batches_scVI(adata, max_considered_genes) -> Tuple[anndata.AnnData, scvi.model.SCVI]:
    """
    Correct batches using scVI. Adds a batch-corrected latent representation to `adata.obsm['X_scVI']`.
    Does internal normalization, so feed raw (filtered) counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data with batch information in `adata.obs['batch']`.
    max_considered_genes : int, optional
        Maximum number of highly variable genes to consider, by default 1000

    Returns
    -------
    adata : anndata.AnnData
        Input data with batch-corrected latent representation in `adata.obsm['X_scVI']`
    model : scvi.model.SCVI
        The trained scVI model
    """
    internal_adata = adata.copy()

    if not hf.is_normalized(internal_adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(internal_adata, target_sum=1e4)
    else:
        vprint("adata is already normalized")

    if max_considered_genes is not "all": # only select HVGs if max_considered_genes is specified as a number
        # do batch aware HVG selection
        print("Selecting highly variable genes...")
        sc.pp.highly_variable_genes(
        internal_adata,
        flavor="seurat_v3",
        n_top_genes=max_considered_genes,      # n_top_genes is the total number of HVGs across all batches
        batch_key="batch"
        )

        n_highly_variable_genes = internal_adata.var['highly_variable'].sum()
        print(f"Found {n_highly_variable_genes} highly variable genes across {internal_adata.obs['batch'].nunique()} batches.")
        
        adata_hvg = internal_adata[:, internal_adata.var['highly_variable']].copy() # only keep HVGs

        #number of genes that are left
        print(f"after applying the boolean mask, there are {adata_hvg.shape[1]} genes left")
        print(f"The shape of adata_hvg is {adata_hvg.shape}")


        # only consider HVGs present in at least 30% of batches (and at very least 2 batches)
        threshold = int(0.3 * internal_adata.obs['batch'].nunique())
        mask = adata_hvg.var.get('highly_variable_nbatches') >= threshold
        adata_hvg = adata_hvg[:, mask].copy()
        print(f"The shape of adata_hvg after thresholding is {adata_hvg.shape}")

        internal_adata = adata_hvg.copy()


    try:
        model = scvi.model.SCVI.load("my_scvi_model/", internal_adata)
    except Exception as e:
        print("No scVI model found. Training a new one...")
        print(f"Model loading failed because of Exception: {e}")
        model = None
        pass

    vprint(model) # check if model is loaded or none

    if model is None:
        # set up model with batch information
        scvi.model.SCVI.setup_anndata(internal_adata, batch_key="batch")

        # Train the model
        model = scvi.model.SCVI(internal_adata)
        model.train()
        model.save("my_scvi_model/", overwrite=True)

    # Get the batch-corrected latent representation (obsm is a matrix like X where each row is a cell and each column is a feature)
    internal_adata.obsm["X_scVI_corrected"] = model.get_normalized_expression()
    print(f"Shape of corrected expression matrix{internal_adata.obsm["X_scVI_corrected"].shape}")

    return internal_adata, model


def correct_batches_scANVI(adata_scvi, pretrained_scVI_model) -> Tuple[anndata.AnnData, scvi.model.SCANVI]:
    """
    Correct batches using scANVI. Adds a batch-corrected latent representation to `adata.obsm['X_scANVI']`.
    Does internal normalization, so feed raw (filtered) counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Input data with batch information in `adata.obs['batch']`.
    pretrained_scVI_model : scvi.model.SCVI
        Pretrained scVI model for transfer learning

    Returns
    -------
    adata : anndata.AnnData
        Input data with batch-corrected latent representation in `adata.obsm['X_scANVI']`
    model : scvi.model.SCANVI
        The trained scANVI model
    """

    internal_adata = adata_scvi.copy()

    if not hf.is_normalized(internal_adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(internal_adata, target_sum=1e4)
    else:
        vprint("adata is already normalized")

    try:
        model = scvi.model.SCANVI.load("my_scanvi_model/", internal_adata)
    except Exception as e:
        print("No scANVI model found. Training a new one...")
        print(f"Model loading failed because of Exception: {e}")
        model = None
        pass

    print(model) # check if model is loaded or none

    if model is None:
        # set up model with batch information
        scvi.model.SCANVI.setup_anndata(internal_adata, batch_key="batch", labels_key="cell_type", unlabeled_category="unlabeled")

        # Train the model
        model = scvi.model.SCANVI.from_scvi_model(scvi_model=pretrained_scVI_model, adata=internal_adata, labels_key="cell_type", unlabeled_category="unlabeled")
        model.train()
        model.save("my_scanvi_model/", overwrite=True)

    # Get the batch-corrected latent representation (obsm is a matrix like X where each row is a cell and each column is a feature)
    internal_adata.obsm["X_scANVI_corrected"] = model.get_normalized_expression()

    return internal_adata, model


def main(input_data_file, output_dir, max_considered_genes):
    adata = sc.read_h5ad(input_data_file)

    adata_scvi, scvi_model = correct_batches_scVI(adata, max_considered_genes=max_considered_genes) 
    corrected_adata, scanvi_model = correct_batches_scANVI(adata_scvi, scvi_model)

    # carry over the corrected latent representations (later get full gene expression and put into layers)
    adata.obsm["X_scVI_corrected"] = corrected_adata.obsm["X_scVI_corrected"]
    adata.obsm["X_scANVI_corrected"] = corrected_adata.obsm["X_scANVI_corrected"]

    # save the processed data to temporary h5ad file, make relevant directory first
    adata.write(os.path.join(output_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_dir, os.path.basename(input_data_file)))



if __name__ == "__main__":

    input_data_file, output_dir, max_considered_genes, verbose = hf.import_cmd_args(4)
    vprint = hf.make_vprint(verbose)

    main(input_data_file, output_dir, max_considered_genes)


