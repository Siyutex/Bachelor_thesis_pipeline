# use this script to gauge suitable filtering thresholds for a given dataset (RUN PER BATCH, ok if you have several similar batches from 1 study)


import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

input_path = r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\OG_data\NCBI\PDAC_non_cancerous\ADJ4"
adata = sc.read_10x_mtx(input_path)

k_MAD = 1 # amount of MADs to add to the cutoffs (3 = industry standard, lower for datasets with very heavy tails)

# violin plot for UMIs per cell
counts_per_cell = np.array(adata.X.sum(axis=1)).flatten() # add all entries for each cell (row)
adata.obs['n_counts'] = counts_per_cell
sc.pl.violin(adata, ['n_counts'], log=True, xlabel="Cells", ylabel="Number of UMI counts per cell")
plt.show()

# violin plot for genes per cell
adata.obs['n_genes'] = adata.X.getnnz(axis=1) # number of entries per cell (including explicitly stored 0s, of which there are none here)
sc.pl.violin(adata, ['n_genes'], log=True, xlabel="Cells", ylabel="Number of genes per cell")
plt.show()

# violin plot for how many cells a gene is expressed in
counts_per_gene = np.array(adata.X.getnnz(axis=0)).flatten() # add all entries for each gene (column)
counts_per_gene = np.sort(counts_per_gene)[::-1]
plt.plot(range(len(counts_per_gene)), counts_per_gene, 'o')
plt.ylabel("Number of cells, where the gene is expressed")
plt.xlabel("Rank of the gene")
plt.title("Distribution of cell counts per gene")
plt.yscale('log')
plt.show()

# violin plot for % mito counts per cell
adata.var["mito"] = adata.var_names.str.startswith("MT-")  # identify mitochondrial genes, assuming they start with "MT-"
adata.obs["pct_counts_mito"] = adata.X[:, adata.var["mito"].values].sum(axis=1) / adata.X.sum(axis=1)
sc.pl.violin(adata, ['pct_counts_mito'], log=False, xlabel="Cells", ylabel="Percentage of mitochondrial counts per cell")
plt.show()


# also show median of each statistic
print(f"Median UMIs per cell: {np.median(counts_per_cell):.1f}")
print(f"Median genes per cell: {np.median(adata.obs['n_genes']):.1f}")
print(f"Median amount of cells a gene is expressed in: {np.median(counts_per_gene):.1f}")
print(f"Median mitochondrial percentage per cell: {np.median(adata.obs['pct_counts_mito']):.3f}")
print("\n")

# also show artihmetic mean of each statistic
print(f"Mean UMIs per cell: {np.mean(counts_per_cell):.1f}")
print(f"Mean genes per cell: {np.mean(adata.obs['n_genes']):.1f}")
print(f"Mean amount of cells a gene is expressed in: {np.mean(counts_per_gene):.1f}")
print(f"Mean mitochondrial percentage per cell: {np.mean(adata.obs['pct_counts_mito']):.3f}")
print("\n")

# suggest thresholds for filtering based on MAD (median absolute deviation)
median_counts_per_cell = np.median(counts_per_cell)
median_genes_per_cell = np.median(adata.obs['n_genes'])
median_cells_per_gene = np.median(counts_per_gene)
median_mito_per_cell = np.median(adata.obs['pct_counts_mito'])

mad_counts_per_cell = np.median(np.abs(counts_per_cell - median_counts_per_cell))
mad_genes_per_cell = np.median(np.abs(adata.obs['n_genes'] - median_genes_per_cell))
mad_cells_per_gene = np.median(np.abs(counts_per_gene - median_cells_per_gene))
mad_mito_per_cell = np.median(np.abs(adata.obs['pct_counts_mito'] - median_mito_per_cell))

print(f"Median absolute deviation (MAD) of UMIs per cell: {mad_counts_per_cell:.1f}")
print(f"Median absolute deviation (MAD) of genes per cell: {mad_genes_per_cell:.1f}")
print(f"Median absolute deviation (MAD) of cells a gene is expressed in: {mad_cells_per_gene:.1f}")
print(f"Median absolute deviation (MAD) of mitochondrial percentage per cell: {mad_mito_per_cell:.3f}")
print("\n")

print("MAD based thresholding:")
print(f"Suggested cutoff for min UMIs per cell: {median_counts_per_cell - 1 * mad_counts_per_cell:.1f}")
print(f"Suggested cutoff for min genes per cell: {median_genes_per_cell - 1 * mad_genes_per_cell:.1f}")
print(f"Suggested cutoff for min cells a gene is expressed in: {median_cells_per_gene - 1 * mad_cells_per_gene:.1f}")
print(f"Suggested cutoff for max mitochondrial percentage per cell: {median_mito_per_cell + 3 * mad_mito_per_cell:.3f}")
print("\n")

print("percentile based thresholding")
print(f"Suggested cutoff for min UMIs per cell: {np.percentile(counts_per_cell, 10):.1f}")
print(f"Suggested cutoff for min genes per cell: {np.percentile(adata.obs['n_genes'], 10):.1f}")
print(f"Suggested cutoff for min cells a gene is expressed in: {np.percentile(counts_per_gene, 10):.1f}")
print(f"Suggested cutoff for max mitochondrial percentage per cell: {np.percentile(adata.obs['pct_counts_mito'], 90):.3f}")
print("\n")
# cells with too many UMIs or genes are going to be removed by doublet detection