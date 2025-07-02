# This script removes trash data from the mtx files, which are used in the biological variance determination pipeline.
# (eg RNAs that are rarely expressed, cells with low quality, etc.)

# the mtx files from NCBI GEO are minimally processed (10x genomics), they contain UMI counts for all cells, including those that are not epithelial cells, and those that are not of high quality. This script removes those cells and RNAs that are not of interest.
# preprocessing similar to paper where data was generated, see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966#:~:text=Summary%20Pancreatic%20ductal%20adenocarcinoma%20,plot%20to%20predict%20the%20overall

import scanpy as sc
import sys
import os

# import mtx file and automatically find tsv files with gene names and cell barcodes
raw_data_dir = sys.argv[1] if len(sys.argv) > 1 else print("Preprocessing: Please provide the path to the raw data directory. It must contain mtx and tsv files from the 10x genomics pipeline") # needed for read_10x_mtx to find the tsv files with gene names and cell barcodes
output_path = sys.argv[2] if len(sys.argv) > 2 else print("Preprocessing: Please provide the path to save the output file.") # path to save the output file, will be saved as h5ad file (AnnData object)

adata = sc.read_10x_mtx(raw_data_dir)  # read mtx file, and tsv files, returns AnnData object


# preprocessing

sc.pp.filter_cells(adata, min_genes=500)  # filter cells with less than 500 genes expressed
sc.pp.filter_genes(adata, min_cells=10)  # filter genes expressed in less than 10 cells

# remove cells with less than 1000 UMI counts
adata.obs['n_counts'] = adata.X.sum(axis=1).A1  # convert sparse matrix to dense array and sum counts per cell
adata = adata[adata.obs['n_counts'] > 1000, :]  # keep cells with more than 1000 UMI counts

# remove cells with high mitochondrial gene expression (>25%)
adata.var["mito"] = adata.var_names.str.startswith("MT-")  # identify mitochondrial genes, assuming they start with "MT-"
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)
# The mitochondrial percentage will be in adata.obs['pct_counts_MT-']
# Filter cells with mitochondrial fraction < 25%
adata = adata[adata.obs['pct_counts_mito'] < 25, :]

# log normalize data
sc.pp.normalize_total(adata, target_sum=1e4)  # normalize each cell to have a total count of 10,000 (eg 1687 UMI counts per cell, actin is 189 UMI counts, so 189/1687 = 0.112, which is 11.2% of the total counts in the cell = 1120 UMI counts per 10,000 UMI counts)
sc.pp.log1p(adata)  # log transform the data

# print to console how many cells and genes are left after preprocessing
print(f"Preprocessing: {adata.shape[0]} cells and {adata.shape[1]} genes left after preprocessing in file: preprocessed_{raw_data_dir} \n")


# save the processed data to h5ad file
os.makedirs(os.path.dirname(output_path), exist_ok=True) # ensure directory exists, if not, create it
adata.write(output_path + ".h5ad")  # save the AnnData object to h5ad file


