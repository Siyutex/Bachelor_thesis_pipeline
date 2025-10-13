import scanpy as sc
import pandas as pd 
import sys











# ------------------------------------------------------------------------------------
adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\batch_corrected\batch_corrected_HVG_PDAC.h5ad")
print(adata.var.columns)

gtf = pd.read_csv(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\gencode.v49.annotation.gtf.gz", sep="\t", comment="#", header=None)

print(gtf.shape)
gtf_genes = gtf[gtf[2] == "gene"]
print(gtf_genes.shape)
print(gtf_genes[2].unique())
gtf_genes = gtf_genes[[0, 3, 4, 8]]
print(gtf_genes.shape)
gtf_genes.columns = ["chromosome", "start", "end", "attributes"]
gtf_genes["gene_id"] = gtf_genes["attributes"].str.extract('gene_id "([^"]+)"')
gtf_genes["gene_name"] = gtf_genes["attributes"].str.extract('gene_name "([^"]+)"')
gtf_genes["gene_type"] = gtf_genes["attributes"].str.extract('gene_type "([^"]+)"')
gtf_genes_protein_coding = gtf_genes[gtf_genes["gene_type"] == "protein_coding"]
gtf_genes_protein_coding = gtf_genes_protein_coding[~gtf_genes_protein_coding["attributes"].str.contains("readthrough")]
gtf_genes_protein_coding = gtf_genes_protein_coding[~gtf_genes_protein_coding["attributes"].str.contains("overlapping_locus")]



# problem: x and y chromosomes contain a lot of genes that are relevant in anndata so don't remove them
# maybe just take the first entry in the assembly for each gene in anndata
# or try to get anndata with ensemnl ids (they are present in the raw data for PDAC and ADJ)
# gene ids do exist but are likely removed in merger in batch correction script

print(gtf_genes_protein_coding.shape)
# save first 100 lines as csv
gtf_small = gtf_genes_protein_coding
gtf_small.to_csv(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\auxiliary_data\annotations\gencode.v49.annotation.gtf.csv", index=False)



# find duplicate labels
import pandas as pd

# Assume gtf_genes_protein_coding and adata are already defined

# ---- 1. Check for duplicates in GTF ----
gtf_dup = gtf_genes_protein_coding["gene_name"].value_counts()
gtf_dup = gtf_dup[gtf_dup > 1]

if not gtf_dup.empty:
    print("Duplicate genes in GTF:")
    for gene, count in gtf_dup.items():
        print(f"  {gene}: {count} entries")

# ---- 2. Check for duplicates in adata ----
adata_dup = adata.var_names.value_counts()
adata_dup = adata_dup[adata_dup > 1]

if not adata_dup.empty:
    print("Duplicate genes in AnnData:")
    for gene, count in adata_dup.items():
        print(f"  {gene}: {count} entries")

# ---- 3. Check for genes in adata that are NOT in GTF ----
gtf_genes_set = set(gtf_genes_protein_coding["gene_name"])
adata_genes_set = set(adata.var_names)
missing_genes = adata_genes_set - gtf_genes_set

if missing_genes:
    print(f"Genes in AnnData not found in GTF ({len(missing_genes)}):")
    print(", ".join(list(missing_genes)[:20]) + ("..." if len(missing_genes) > 20 else ""))

sys.exit()

coords_matched = gtf_genes.set_index("gene_name").reindex(adata.var_names)

"""adata.var["chromosome"] = coords_matched["chromosome"].values
adata.var["start"] = coords_matched["start"].values
adata.var["end"] = coords_matched["end"].values

print(adata.var_names[:10])
print(adata.var["chromosome"][:10])"""