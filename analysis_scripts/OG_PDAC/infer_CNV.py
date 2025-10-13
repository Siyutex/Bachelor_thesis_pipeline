# infer copy number variations on an aggregated / batch corrected h5ad file
# Input is an aggregated (batch corrected) h5ad file
# output is an h5ad, just containing ductal cells, with the following new annotations:
# adata.layers["gene_values_cnv"] = number of copies per gene
# adata.obs["cnvs"] = summed number of gene copies for all genes per cell (cancer cells should be more than normal cells)
# adata.obsm["X_cnv"] = smoothed and denoised gene expression along genomic locations

import infercnvpy as cnv
import helper_functions as hf
import scanpy as sc
import os
import pandas as pd
import sys
import numpy as np
from scipy import sparse

if __name__ == "__main__":
    # import cmd args
    input_data_file, output_data_dir, refernce_genome_path, corrected_representation, verbose = hf.import_cmd_args(5)
    vprint = hf.make_vprint(verbose)

    PREFIX = "batch_corrected_HVG_"

    # import adata
    adata = sc.read_h5ad(input_data_file)
    # use batch corrected representation
    if corrected_representation is not None:
        adata.X = adata.obsm[corrected_representation].copy()

    # isolate ductal cells to save time
    vprint("Isolating ductal cells...")
    ductal_cells = adata.obs["cell_type"] == "ductal_cell"
    adata = adata[ductal_cells, :]

    # preprocess (normalization and logarithmization expected by inferCNVpy)
    if not hf.is_normalized(adata):
        vprint("Normalizing adata...")
        sc.pp.normalize_total(adata, target_sum=1e4)
    else:
        vprint("adata is already normalized")
    vprint("Logarithmizing data...")
    sc.pp.log1p(adata)

    # annotate gene coordinates
    vprint("Annotating gene coordinates...")
    ref_genome = pd.read_csv(refernce_genome_path, sep="\t", comment="#", header=None)
    ref_genome = ref_genome[ref_genome[2] == "gene"] # only interested in genes (exons, transcripts etc not relevant for CNV inference)
    ref_genome = ref_genome.iloc[:, [0, 3, 4, 8]] # limit to columns of interest
    ref_genome.columns = ["chromosome", "start", "end", "attributes"] # column names matched according to documention of reference genome
    ref_genome["gene_id"] = ref_genome["attributes"].str.extract('gene_id "([^"]+)"') # add new column with ensembl id
    ref_genome["gene_id"] = ref_genome["gene_id"].str.replace(r"\.\d+$", "", regex=True) # Remove the period and digits at the end of each string (this is a version number, it is often removed in public data, so we remove it here to avoid conflicts)

    print("Matching index of adata to reference genome...")
    coords_matched = ref_genome.set_index("gene_id").reindex(adata.var_names)

    print("Adding chromosomal coordinates to adata...")
    adata.var["chromosome"] = coords_matched["chromosome"].values
    adata.var["start"] = coords_matched["start"].values
    adata.var["end"] = coords_matched["end"].values

    # infer CNV
    vprint("Inferring CNV...")
    cnv.tl.infercnv(adata, reference_key="cancer_state", reference_cat="non_cancerous", calculate_gene_values=True, key_added="cnv") # reference of which cells are normal and should not have CNVs (from annotated input data)

    # add CNVs per cell to obs
    vprint("Adding CNVs per cell to obs...")
    X = adata.layers["gene_values_cnv"]

    if sparse.issparse(X):
        summed = np.array(X.sum(axis=1)).ravel()   # convert from matrix to 1D ndarray
    else:
        summed = np.nansum(X, axis=1)              # ignore NaNs safely

    adata.obs["summed_cnvs"] = summed

    print(adata.obs.keys())
    print(adata.obs["summed_cnvs"][:5])

    # save results
    vprint("Saving results...")
    adata.write(os.path.join(output_data_dir, f"cnv_{os.path.basename(input_data_file).removeprefix(PREFIX)}"), compression="gzip")
    print("Output: " + os.path.join(output_data_dir, f"cnv_{os.path.basename(input_data_file).removeprefix(PREFIX)}"))
