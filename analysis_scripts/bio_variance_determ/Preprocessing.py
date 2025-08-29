# removes:
"""
Cells with less than 500 genes expressed
Genes expressed in less than 10 cells
Cells with less than 1000 UMI counts
Cells with high mitochondrial gene expression (>25%)
Genes not calssified as "highly variable" 
"""
# and log normalizes the data (per 10K counts, log1p)

# the mtx files from NCBI GEO are minimally processed (10x genomics), they contain UMI counts for all cells, including those that are not epithelial cells, and those that are not of high quality. This script removes those cells and RNAs that are not of interest.
# preprocessing similar to paper where data was generated, see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966#:~:text=Summary%20Pancreatic%20ductal%20adenocarcinoma%20,plot%20to%20predict%20the%20overall

import scanpy as sc
import sys
import os
import tempfile
import tarfile
import pandas as pd

# import command line arguments from ececutor script
input_data_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the raw data directory. It must contain mtx and tsv files from the 10x genomics pipeline") # can be a folder or a file, depending on the datatype
output_data_path = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to where the output should be saved")
input_data_type = sys.argv[3] if len(sys.argv) > 3 else print("Please provide the data type. It must be one of the following: MTX_TSVs_in_subfolders, compressed_MTX_TSVs, dot_matrix_files") # string datatype variable to indicate which pipeline mode was chosen in the executor script


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
        adata = sc.read_10x_mtx(temp_dir)
elif input_data_type == "dot_matrix_files":
    # directly read the forwarded .matrix file = raw_data
    temp_df = pd.read_csv(input_data_path, sep="\t", index_col=0) # rows are genes, columns are cells, need to transpose to fit with annData format
    adata = sc.AnnData(temp_df.T) # transpose the dataframe to have cells as rows and genes as columns



# preprocessing

sc.pp.filter_cells(adata, min_genes=500)  # filter cells with less than 500 genes expressed
sc.pp.filter_genes(adata, min_cells=int(0.01*adata.n_obs))  # filter genes expressed in less than 1% of cells

# remove cells with less than 1000 UMI counts
adata.obs['n_counts'] = adata.X.sum(axis=1)  # convert sparse matrix to dense array and sum counts per cell (axis=1 means it looks through all columns(genes) per row(cell))
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

'''# isolate highly variable genes (harmony batch correction works in PCA space, meaning non HVGs are not considered)
# also for pseudotime and the final GRN, only considering HVGs is beneficial so we remove all other right here
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var['highly_variable']]'''

# print to console how many cells and genes are left after preprocessing
print(f"{adata.shape[0]} cells and {adata.shape[1]} genes left after preprocessing the file: {os.path.basename(input_data_path)}")


# save the processed data to temporary h5ad file, make relevant directory first
os.makedirs(os.path.join(output_data_path, "preprocessing"),exist_ok=True)
final_output_path = os.path.join(output_data_path, "preprocessing", os.path.basename(input_data_path) + ".h5ad")
adata.write(final_output_path)

# foward the path to the temporary file to the executor script via stdout
print("Output: " + final_output_path)


