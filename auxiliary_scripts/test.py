import helper_functions as hf
vprint = hf.make_vprint(True)
import numpy as np
import scanpy as sc


adata = sc.read_h5ad(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\test_CNV_PDAC_ductal_cell.h5ad")

def add_additional_annotations(adata):
    # set adata.obs["cancer_state_inferred"] to cancerous or non_cancerous, based on cnv score percentile
    # define percentile based on expected amount of cancerous cells from cancer_state (batch origin)
    percentile = (len(adata[adata.obs["cancer_state"] == "non_cancerous"]) / adata.shape[0]) * 100 #percentile is a float between 0 and 100
    adata.obs["cancer_state_inferred"] = np.where(adata.obs["cnv_score"] > np.percentile(adata.obs["cnv_score"], percentile), "cancerous", "non_cancerous")


add_additional_annotations(adata)
print(adata)

adata.write(r"C:\Users\Julian\Documents\not_synced\Github\Bachelor_thesis_pipeline\Data\output_storage\CNV\test_CNV_PDAC_ductal_cell.h5ad", compression="gzip")