import scanpy as sc
import sys
import numpy as np
import pandas as pd

# import gz compressed h5ad file (AnnData object)
ADJ_paths = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the gz compressed h5ad file.")  # sys.argv takes command line arguments (or arguments from other scripts, here: biological_variance_pipeline_executor.py), first is the script name, second is the input file path
PDAC_paths = sys.argv[2] if len(sys.argv) > 2 else print("Please provide the path to save the output file.")  # path to save the output file, will be saved as gz compressed h5ad file (AnnData object)

ADJ_paths = ADJ_paths.split(",")  # split the input paths by comma to get a list of paths
PDAC_paths = PDAC_paths.split(",")



def find_common_gene_set():
    all_paths = ADJ_paths + PDAC_paths  # Combine both paths into a list
    common_genes = None  # Initialize common genes variable
    for f in all_paths:
        print(f"Variance: Processing file {f} to find common genes...")
        adata = sc.read_h5ad(f)  # Read each h5ad file
        if common_genes is None:
            common_genes = set(adata.var_names)  # Initialize with the first file's genes
        else:
            common_genes.intersection_update(adata.var_names)  # Find intersection with subsequent files
    print(f"Variance: Common genes found: {len(common_genes)}")
    return common_genes  # Return the set of common genes across all files


def create_mean_expr_DF(file_paths, common_genes):

    # Container for pseudobulk data
    pseudobulk_means = []


    for f in file_paths:
        adata = sc.read_h5ad(f)
        adata = adata[:, adata.var_names.isin(common_genes)]  # Filter to keep only common genes
        
        # Extract patient ID (optional: parse from filename)
        patient_id = f.split("/")[-1].replace(".h5ad", "") # split path into strings with "/" as delimiter, take last string (file name), remove ".h5ad" suffix
        
        expr = adata.X # Extract expression matrix from AnnData object
            
        # If expr is sparse matrix, convert to dense
        if not isinstance(expr, np.ndarray):
            expr = expr.toarray()
        
        # Mean expression per gene across cells
        gene_means = expr.mean(axis=0)
        
        # Store as DataFrame row
        pseudobulk_means.append(pd.Series(gene_means, index=adata.var_names, name=patient_id))

    # Combine into final matrix: patients x genes
    pseudobulk_df = pd.DataFrame(pseudobulk_means)
    print(pseudobulk_df.shape)  # (6, num_genes)
    return pseudobulk_df  # Return the DataFrame containing mean expression values



def compute_variance(pseudobulk_df):
    """Compute variance between all genes across all patients."""

    # Compute variance across patients (rows) for each gene (columns)
    variance_df = pseudobulk_df.var(axis=0)
    stdev_df = pseudobulk_df.std(axis=0)
    
    # Sort genes by variance (descending)
    sorted_variance = variance_df.sort_values(ascending=False)
    sorted_stdev = stdev_df.sort_values(ascending=False)

    # compute average variance across all genes and print it
    avg_variance = sorted_variance.mean()
    avg_stdev = sorted_stdev.mean()
    print(f"Variance: Average variance across all genes: {avg_variance}")
    print(f"Variance: Standard deviation across all genes: {avg_stdev} \n")



# execute the functions
common_genes = find_common_gene_set()  # Find common gene set across all files


# run variance computation for both ADJ and PDAC files
ADJ_DF = create_mean_expr_DF(ADJ_paths, common_genes)
PDAC_DF = create_mean_expr_DF(PDAC_paths, common_genes)
compute_variance(ADJ_DF)
compute_variance(PDAC_DF)