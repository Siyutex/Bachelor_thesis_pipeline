# infer copy number variations on an aggregated / batch corrected h5ad file
# Input is an aggregated (batch corrected) h5ad file
# output is an h5ad, just containing ductal cells, with the following new annotations:
# adata.obsm["{representation}_gene_values_cnv"] = number of copies per gene
# adata.obs["summed_cnvs"] = summed number of gene copies for all genes per cell (cancer cells should be more than normal cells)
# adata.obsm["{representation}_cnv"] = smoothed and denoised gene expression along genomic locations
# adata.obs["cnv_score"] = cnv score of leiden clusters (see inferCNVpy.tl.cnv_score)

import infercnvpy as cnv
import helper_functions as hf
import scanpy as sc
import os
import pandas as pd
import numpy as np
from scipy import sparse

def check_assign_corrected_representation(adata, corrected_representation):

    # use batch corrected representation if given, directly modify passed adata
    if corrected_representation is not None and corrected_representation in adata.obsm.keys():
        vprint(f"Using corrected representation {corrected_representation}...")
        adata.X = adata.obsm[corrected_representation].copy()
    elif corrected_representation is not None and corrected_representation not in adata.obsm.keys():
        raise ValueError(f"Corrected representation {corrected_representation} not found in adata.obsm")
    elif corrected_representation is None:
        vprint("Using adata.X...")


def check_assign_cell_type(adata, cell_type) -> sc.AnnData:

    # limit adata to chosen cell_type if given
    if cell_type is not None and cell_type in adata.obs["cell_type"].values:
        vprint(f"Isolating {cell_type}...")
        adata = adata[adata.obs["cell_type"] == cell_type, :].copy()
    elif cell_type is not None and cell_type not in adata.obs["cell_type"].values:
        raise ValueError(f"Cell type {cell_type} not found in adata.obs['cell_type']")
    elif cell_type is None:
        vprint("Using all cells...") 

    return adata # we need to return the object here because slicing creates a new object, so the adata passed as a paremeter is not modified


def preprocess_data(adata):

    # preprocess (normalization and logarithmization expected by inferCNVpy)
    if not hf.is_normalized(adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(adata, target_sum=1e4)
    else:
        vprint("adata is already normalized")

    vprint("Logarithmizing data...") # minimally, scuffed bcs no check, but in this pipeline there is no log1p before this and we only ever tranform X on internal adata
    sc.pp.log1p(adata)


def assign_chromosomal_coordinats(adata, reference_genome_path):

    # annotate gene coordinates
    vprint("Preparing reference genome...")
    ref_genome = pd.read_csv(refernce_genome_path, sep="\t", comment="#", header=None)
    ref_genome = ref_genome[ref_genome[2] == "gene"] # only interested in genes (exons, transcripts etc not relevant for CNV inference)
    ref_genome = ref_genome.iloc[:, [0, 3, 4, 8]] # limit to columns of interest
    ref_genome.columns = ["chromosome", "start", "end", "attributes"] # column names matched according to documention of reference genome
    ref_genome["gene_id"] = ref_genome["attributes"].str.extract('gene_id "([^"]+)"') # add new column with ensembl id
    ref_genome["gene_id"] = ref_genome["gene_id"].str.replace(r"\.\d+$", "", regex=True) # Remove the period and digits at the end of each string (this is a version number, it is often removed in public data, so we remove it here to avoid conflicts)

    vprint("Matching index of adata to reference genome...")
    coords_matched = ref_genome.set_index("gene_id").reindex(adata.var_names)

    vprint("Adding chromosomal coordinates to adata...")
    adata.var["chromosome"] = coords_matched["chromosome"].values
    adata.var["start"] = coords_matched["start"].values
    adata.var["end"] = coords_matched["end"].values

def add_additional_annotations(adata):

    # add CNVs per cell to obs
    vprint("Adding CNVs per cell to adata.obs...")
    X = adata.layers["gene_values_cnv"] # gene values cnv will have many nans if the used geneset does not densly cover all chromosomes (sliding window must have at least n genes per window, or it skips that window)

    if sparse.issparse(X):
        summed = np.nansum(np.array(X), axis=1) # convert from matrix to 1D ndarray
    else:
        summed = np.nansum(X, axis=1)              # ignore NaNs safely

    adata.obs["summed_cnvs"] = summed

def main(input_data_file, output_data_dir, refernce_genome_path, corrected_representation, cell_type):

    # import adata
    print("Reading data...")
    adata = sc.read_h5ad(input_data_file)
    internal_adata = adata.copy()

    # check if corrected representation should be used, if so modify adata
    print("Checking if corrected representation should be used...")
    check_assign_corrected_representation(internal_adata, corrected_representation)

    # check if cell type should be limited, if so modify adata
    print("Checking if cell type should be limited...")
    internal_adata = check_assign_cell_type(internal_adata, cell_type)
    adata = check_assign_cell_type(adata, cell_type) # we have also limit the actual adata to the cell type, if it is passed so the obsm annotation fits

    vprint(f"Internal adata size: {internal_adata.shape}")
    vprint(f"Original adata size: {adata.shape}")

    # preprocess (normalization and log1p required by inferCNVpy)
    print("Preprocessing data...")
    preprocess_data(internal_adata)

    # assign chromosomal coordinates to genes
    print("Assigning chromosomal coordinates...")
    assign_chromosomal_coordinats(internal_adata, refernce_genome_path)

    # infer CNVs
    print("Inferring CNV...")
    cnv.tl.infercnv(internal_adata, reference_key="cancer_state", reference_cat="non_cancerous", calculate_gene_values=True, key_added="cnv", chunksize=100) # reference of which cells are normal and should not have CNVs (from annotated input data)

    # get cnv score for each cluster of cells (high score = prbly aneuploidy)
    # first cnv: pca, neighbours, leiden (all of these use scanpy default values in the underlying function)
    cnv.tl.pca(internal_adata) # uses obsm["X_{use_rep}"] -> so just the smoothed, denoised, cell comparable log fold changes in gene expression over the reference
    cnv.pp.neighbors(internal_adata)
    cnv.tl.leiden(internal_adata)

    cnv.tl.cnv_score(internal_adata) # mean absolute values of X_cnv of each cluster (so mean across all cells from that cluster)

    # add additional annotations
    print("Adding additional annotations...")
    add_additional_annotations(internal_adata)

    # assign new annotations to original adata
    if corrected_representation is not None:
        adata.obsm[f"{corrected_representation}_cnv"] = internal_adata.obsm["X_cnv"]
        adata.obsm[f"{corrected_representation}_gene_values_cnv"] = internal_adata.layers["gene_values_cnv"]
        adata.obs["summed_cnvs"] = internal_adata.obs["summed_cnvs"]
        adata.obs["cnv_score"] = internal_adata.obs["cnv_score"]
    else:
        adata.obsm["X_cnv"] = internal_adata.obsm["X_cnv"]
        adata.obsm["X_gene_values_cnv"] = internal_adata.layers["gene_values_cnv"] # still add to obsm so downstream processing is uniform
        adata.obs["summed_cnvs"] = internal_adata.obs["summed_cnvs"]
        adata.obs["cnv_score"] = internal_adata.obs["cnv_score"]

    # save results
    print("Saving results...")
    adata.write(os.path.join(output_data_dir, os.path.basename(input_data_file)), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, os.path.basename(input_data_file)))






if __name__ == "__main__":
    # import cmd args
    input_data_file, output_data_dir, refernce_genome_path, corrected_representation, cell_type, verbose = hf.import_cmd_args(6)
    vprint = hf.make_vprint(verbose)

    main(input_data_file, output_data_dir, refernce_genome_path, corrected_representation, cell_type)




    


    

