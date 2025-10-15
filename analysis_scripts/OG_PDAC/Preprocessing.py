# Filtering and doublet removal, default thresholds visible in filtering_params dataclass

# log normalization and HVG selection should be done later when needed in other scripts
# (scVI needs raw counts, so normalization here would brick that model)
# so the output is raw, filtered counts

import scanpy as sc
import os
import tempfile
import tarfile
import pandas as pd
import numpy as np
from scipy import sparse
import helper_functions as hf


def read_input_data(input_data_type, input_data_path, var_names):
    """
    Read different types of single-cell input data formats into an AnnData object.
    """
    if input_data_type == "MTX_TSVs_in_subfolders":
        # Directly read 10X-formatted data (uncompressed)
        return sc.read_10x_mtx(input_data_path, var_names=var_names)

    elif input_data_type == "compressed_MTX_TSVs_in_subfolders":
        return _read_compressed_10x(input_data_path, var_names)

    elif input_data_type == "dot_matrix_files":
        return _read_dot_matrix(input_data_path, var_names)

    else:
        raise ValueError(f"Unsupported input_data_type: {input_data_type}")


def _read_compressed_10x(input_dir, var_names):
    """
    Extract .tar.gz archives to a temporary directory and load them using scanpy.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        for file in os.listdir(input_dir):
            if file.endswith((".tar.gz", ".tgz")):
                file_path = os.path.join(input_dir, file)
                with tarfile.open(file_path, "r:gz") as tar:
                    tar.extractall(path=temp_dir)
        extracted_dirs = [os.path.join(temp_dir, d) for d in os.listdir(temp_dir)]
        if not extracted_dirs:
            raise FileNotFoundError(f"No valid extracted folder found in {input_dir}")
        return sc.read_10x_mtx(extracted_dirs[0], var_names=var_names)


def _read_dot_matrix(file_path, var_names):
    """
    Read a .matrix file (genes x cells) and convert to AnnData.
    """
    if var_names != "gene_symbols":
        raise ValueError(
            "var_names must be 'gene_symbols' when input_data_type is 'dot_matrix_files'. Set use_ensembl_ids to False in pipeline_executor."
        )

    df = pd.read_csv(file_path, sep="\t", index_col=0)
    return sc.AnnData(df.T)


def filter_cells_genes(adata, min_n_genes_percentile, min_n_cells_percentage):

    # standardize adata to use scipy sparse matrix for X, or getnnz will not work
    if type(adata.X) != sparse.csc_matrix:
        adata.X = sparse.csc_matrix(adata.X)

    # filtering cells and genes
    vprint("filtering cells with low gene counts")
    adata.obs['n_genes'] = adata.X.getnnz(axis=1) # number of entries per cell (including explicitly stored 0s, of which there are none here)
    sc.pp.filter_cells(adata, min_genes=int(np.percentile(adata.obs['n_genes'], min_n_genes_percentile))) # filter cells that have less genes expressed than the 10th percentile (the cell with the highest amount of genes in the bottom 10% of cells)
    
    vprint("filtering genes with low cell counts")
    sc.pp.filter_genes(adata, min_cells=int(min_n_cells_percentage*adata.n_obs))  # filter genes expressed in less than 1% of cells


def filter_UMI_counts(adata, min_n_UMIs_percentile):

    # remove cells with less UMI counts than the 10th percentile
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten() # add all umi counts for each cell (row)
    adata = adata[adata.obs['n_counts'] > int(np.percentile(adata.obs['n_counts'], min_n_UMIs_percentile)), :]  # keep cells with more UMI counts than the 10th percentile


def filter_mito_percentage(adata, max_n_MADs, var_names):

    if var_names == "gene_ids":
        adata.var["mito"] = adata.var["gene_symbols"].str.startswith("MT-")  # identify mitochondrial genes, assuming they start with "MT-"
    elif var_names == "gene_symbols":
        adata.var["mito"] = adata.var_names.str.startswith("MT-")
    
    adata.obs["pct_counts_mito"] = adata.X[:, adata.var["mito"].values].sum(axis=1) / adata.X.sum(axis=1)
    mito_cutoff = np.median(adata.obs['pct_counts_mito']) + max_n_MADs * np.median(np.abs(adata.obs['pct_counts_mito'] - np.median(adata.obs['pct_counts_mito']))) # median + MAD
    adata = adata[adata.obs['pct_counts_mito'] < mito_cutoff, :]


def filter_doublets(adata, expected_doublet_percentage):

    # remove doublets using scrublet
    sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_percentage) # boolean prediction in .obs['predicted_doublet']
    adata = adata[~adata.obs['predicted_doublet']] # ~ is a bitwise NOT operator, so we keep all cells where predicted_doublet == False


def save_output(adata, input_data_file_or_dir, output_dir):
    # save the processed data to temporary h5ad file, make relevant directory first
    final_output_path = os.path.join(output_dir, os.path.basename(input_data_file_or_dir) + ".h5ad")
    adata.write(final_output_path, compression="gzip")

    # foward the path to the temporary file to the executor script via stdout
    print("Output: " + final_output_path)


def main(input_data_file_or_dir, 
         output_dir, 
         input_data_type, 
         var_names, 
         min_n_genes_percentile, 
         min_n_cells_percentage, 
         min_n_UMIs_percentile, 
         max_n_MADs, 
         expected_doublet_percentage):

    print("reading input data...")
    adata = read_input_data(input_data_type, input_data_file_or_dir, var_names)

    print(f"{adata.shape[0]} cells and {adata.shape[1]} genes present before filtering")
    print("filtering cells and genes...")
    filter_cells_genes(adata, min_n_genes_percentile, min_n_cells_percentage)

    print("filtering cells with low UMI counts...")
    filter_UMI_counts(adata, min_n_UMIs_percentile)

    print("filtering cells with high mitochondrial gene expression...")
    filter_mito_percentage(adata, max_n_MADs, var_names=var_names)

    print("removing doublets...")
    filter_doublets(adata, expected_doublet_percentage)

    # print to console how many cells and genes are left after preprocessing
    print(f"{adata.shape[0]} cells and {adata.shape[1]} genes left after preprocessing")
    print(f"Anndata var columns: {adata.var.columns}")

    print("Saving output...")
    save_output(adata, input_data_file_or_dir, output_dir)










if __name__ == "__main__":

    ## get input arguments
    print("Importing command line arguments...")
    input_data_file_or_dir, output_dir, input_data_type, filtering_params_list, use_ensembl_ids, verbose = hf.import_cmd_args(6)
    vprint = hf.make_vprint(verbose)

    # set var_names
    if use_ensembl_ids:
        var_names = "gene_ids"
    else:
        var_names = "gene_symbols"

    # set filtering params
    min_n_genes_percentile = filtering_params_list[0]
    min_n_cells_percentage = filtering_params_list[1]
    min_n_UMIs_percentile = filtering_params_list[2]
    max_n_MADs = filtering_params_list[3]
    expected_doublet_percentage = filtering_params_list[4]


    ## execute functions
    main(input_data_file_or_dir, 
         output_dir, input_data_type, 
         var_names,
         min_n_genes_percentile, 
         min_n_cells_percentage, 
         min_n_UMIs_percentile, 
         max_n_MADs, 
         expected_doublet_percentage)

