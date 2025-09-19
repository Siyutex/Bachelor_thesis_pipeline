# this scrript checks if a UMI matrix is raw features or filtered features 
# eg GDC gives raw UMI counts, CancerSCEM gives filtered UMI counts (which means, not all barcodes are included, like empty droplets)

import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tempfile
import os
import tarfile

def check_matrix(filepath, type):

    if type == "CancerSCEM":
        temp_df = pd.read_csv(filepath, sep="\t", index_col=0) # rows are genes, columns are cells, need to transpose to fit with annData format
        adata = sc.AnnData(temp_df.T) # transpose the dataframe to have cells as rows and genes as columns
    if type == "GDC":
        with tempfile.TemporaryDirectory() as temp_dir:
            for file in os.listdir(filepath):
                if file.endswith(".tar.gz") or file.endswith(".tgz"):
                    file_path = os.path.join(filepath, file)
                    with tarfile.open(file_path, "r:gz") as tar:
                        tar.extractall(path=temp_dir)
            adata = sc.read_10x_mtx(os.path.join(temp_dir, os.listdir(temp_dir)[0]))
    if type == "NCBI":
        # read mtx file, and tsv files from the current folder in the raw_data directory
        adata = sc.read_10x_mtx(filepath)  
    
    
    print(f" Shape of adata is {adata.shape}")
    
    # Compute total counts per barcode (cell)
    counts_per_cell = np.array(adata.X.sum(axis=1)).flatten()
    n_barcodes = len(counts_per_cell)
    
    print(f"✅ Loaded {filepath}")
    print(f"   Number of barcodes: {n_barcodes}")
    print(f"   Median UMIs per barcode: {np.median(counts_per_cell):.1f}")
    print(f"   Max UMIs per barcode: {np.max(counts_per_cell):.0f}")
    
    # Sort and plot knee plot
    counts_sorted = np.sort(counts_per_cell)[::-1]
    plt.figure(figsize=(6,4))
    plt.plot(counts_sorted)
    plt.yscale("log")
    plt.xlabel("Barcodes sorted by total UMIs")
    plt.ylabel("Total UMIs per barcode (log scale)")
    plt.title("UMI Knee Plot")
    plt.show()
    
    # Heuristic: decide raw vs filtered
    if n_barcodes > 20000:
        print("⚠️ Likely RAW matrix (includes empty droplets).")
    elif n_barcodes < 15000 and np.median(counts_per_cell) > 500:
        print("ℹ️ Likely FILTERED matrix (cell-called).")
    else:
        print("❓ Ambiguous — check the knee plot.")




for subdir in os.listdir(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\OG_data\NCBI\PDAC_cancerous"):
    filepath = os.path.join(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\OG_data\NCBI\PDAC_cancerous", subdir)
    check_matrix(filepath, "NCBI")
