# We assume the input is the OG PDAC dataset. It was sequenced with 10x chromium, 3000 recovered cells per batch
# thus we expect ~2.4% doublets according to https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
# removes:
"""
Cells with less than 500 genes expressed
Genes expressed in less than 1% of cells
Cells with less than 1000 UMI counts
Cells with high mitochondrial gene expression (>25%) 
correct for doublets
"""
# log normalization and HVG selection should be done later when needed in other scripts
# (scVI needs raw counts, so normalization here would brick that model)
# so the output is raw, filtered counts

# the mtx files from NCBI GEO are minimally processed (10x genomics), they contain UMI counts for all cells, including those that are not epithelial cells, and those that are not of high quality. This script removes those cells and RNAs that are not of interest.
# preprocessing similar to paper where data was generated, see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966#:~:text=Summary%20Pancreatic%20ductal%20adenocarcinoma%20,plot%20to%20predict%20the%20overall

import scanpy as sc
import sys
import os
import tempfile
import tarfile # needed to extract compressed tar files
import pandas as pd
import numpy as np


# import command line arguments from ececutor script
input_data_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the raw data directory. It must contain mtx and tsv files from the 10x genomics pipeline") # can be a folder or a file, depending on the datatype
output_data_path = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to where the output should be saved")
input_data_type = sys.argv[3] if len(sys.argv) > 3 else print("Please provide the data type. It must be one of the following: MTX_TSVs_in_subfolders, compressed_MTX_TSVs, dot_matrix_files") # string datatype variable to indicate which pipeline mode was chosen in the executor script


# read data into adata, depending on input data type
print("Reading data")
if input_data_type == "MTX_TSVs_in_subfolders":
    # read mtx file, and tsv files from the current folder in the raw_data directory
    adata = sc.read_10x_mtx(input_data_path)  
elif input_data_type == "compressed_MTX_TSVs_in_subfolders":
    # extract compressed files to a temporary directory and read them from there
    with tempfile.TemporaryDirectory() as temp_dir:
        for file in os.listdir(input_data_path):
            if file.endswith(".tar.gz") or file.endswith(".tgz"):
                file_path = os.path.join(input_data_path, file)
                with tarfile.open(file_path, "r:gz") as tar:
                    tar.extractall(path=temp_dir)
        adata = sc.read_10x_mtx(os.path.join(temp_dir, os.listdir(temp_dir)[0]))
elif input_data_type == "dot_matrix_files":
    # directly read the forwarded .matrix file = raw_data
    temp_df = pd.read_csv(input_data_path, sep="\t", index_col=0) # rows are genes, columns are cells, need to transpose to fit with annData format
    adata = sc.AnnData(temp_df.T) # transpose the dataframe to have cells as rows and genes as columns

print(f"{adata.shape[0]} cells and {adata.shape[1]} genes present before preprocessing")

print("filtering cells with low gene counts")
print("filtering genes with low cell counts")
# filtering cells and genes
adata.obs['n_genes'] = adata.X.getnnz(axis=1) # number of entries per cell (including explicitly stored 0s, of which there are none here)
sc.pp.filter_cells(adata, min_genes=int(np.percentile(adata.obs['n_genes'], 10))) # filter cells that have less genes expressed than the 10th percentile (the cell with the highest amount of genes in the bottom 10% of cells)
sc.pp.filter_genes(adata, min_cells=int(0.01*adata.n_obs))  # filter genes expressed in less than 1% of cells

print("removing low UMI count cells")
# remove cells with less UMI counts than the 10th percentile
adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten() # add all umi counts for each cell (row)
adata = adata[adata.obs['n_counts'] > int(np.percentile(adata.obs['n_counts'], 10)), :]  # keep cells with more UMI counts than the 10th percentile

print("removing cells with high mitochondrial gene expression")
adata.var["mito"] = adata.var_names.str.startswith("MT-")  # identify mitochondrial genes, assuming they start with "MT-"
adata.obs["pct_counts_mito"] = adata.X[:, adata.var["mito"].values].sum(axis=1) / adata.X.sum(axis=1)
mito_cutoff = np.median(adata.obs['pct_counts_mito']) + np.median(np.abs(adata.obs['pct_counts_mito'] - np.median(adata.obs['pct_counts_mito']))) # median + MAD
adata = adata[adata.obs['pct_counts_mito'] < mito_cutoff, :]

# remove doublets 
print("removing doublets")
sc.pp.scrublet(adata, expected_doublet_rate=0.024) # boolean prediction in .obs['predicted_doublet']
adata = adata[~adata.obs['predicted_doublet']] # ~ is a bitwise NOT operator, so we keep all cells where predicted_doublet == False

# print to console how many cells and genes are left after preprocessing
print(f"{adata.shape[0]} cells and {adata.shape[1]} genes left after preprocessing")


# save the processed data to temporary h5ad file, make relevant directory first
final_output_path = os.path.join(output_data_path, os.path.basename(input_data_path) + ".h5ad")
adata.write(final_output_path, compression="gzip")

# foward the path to the temporary file to the executor script via stdout
print("Output: " + final_output_path)


