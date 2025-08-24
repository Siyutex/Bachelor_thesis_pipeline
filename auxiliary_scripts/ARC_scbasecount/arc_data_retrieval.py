# Metadata file: "gs://arc-scbasecount/2025-02-25/metadata/GeneFull/Homo_sapiens/sample_metadata.parquet"

# 16077 unique accession numbers, no duplicates

# relevant valid disease types (only for samples without perturbation (according to no perturbation list in function above)): 
relevant_valid_disease_types = [
    "gliomablastoma", "glioblastoma (GBM)",
    "prostate cancer", "prostate adenocarcinoma",
    "lung cancer", "lung adenocarcinoma",
    "chronic myelomonocytic leukemia", "chronic myelomonocytic leukemia (CMML)", "acute lymphoblastic leukemia",
    "pancreatic ductal adenocarcinoma", "pancreatic ductal adenocarcinoma (PDAC)",
    "colorectal cancer",
    "breast cancer",
    "gastric cancer",
    "testicular cancer",
    "bladder cancer"
]

# relevant perturbation types that count as "no perturbation" according to personal inspection
no_perturbation = ["none", "none reported", "none specified", "none specified (healthy)", "none specified; cohort labeled as 'Control'", "none (control group)", "none specified; study focuses on genetic variations", "none (control biopsy)", "none specified (healthy control)", "none specified (control group)", "none specified; samples obtained from patients undergoing surgical resection"]


import pandas as pd
import os

## INITIALIZE THE SCRIPT

# relevant directories for this script
script_dir = os.path.dirname(__file__)
data_dir = os.path.join(script_dir, "..", "..", "auxiliary_data", "ARC_scbasecount")
# load the data into a pandas dataframe and print names of all columns
df = pd.read_parquet(os.path.join(data_dir, "GeneFull_HomoSapiens_sample_metadata.parquet"), engine="pyarrow")



## FUNCTIONS

#get all possible values for the tissue and disease parameter and save in a list (for the samples without perturbation)
def get_valid_paramters(all_perturbations: bool=False):

    tissues = []
    diseases = []
    perturbations = []

    for row in range(df.shape[0]):

        # "none" is a valid perturbation value (I tested this seperately)
        if df.loc[row, "perturbation"] in no_perturbation or all_perturbations == True:
            if df.loc[row, "tissue"] not in tissues:
                tissues.append(df.loc[row, "tissue"])
            if df.loc[row, "disease"] not in diseases:
                diseases.append(df.loc[row, "disease"])
            if df.loc[row, "perturbation"] not in perturbations:
                perturbations.append(df.loc[row, "perturbation"])

    return tissues, diseases, perturbations

# used this to search for disease valid and relevant disease types 
'''diseases = get_valid_paramters()[1]
searchterm = "glio"

counter = 0
for element in diseases:
    if searchterm in element:
        print(element + "\n")
        counter += 1
print(f"Number of disease types containing the searchterm '{searchterm}': {counter}")'''

# create a new dataframe that only contains rows that have certain disease type (cancer) and perturbation type (none)
def reduce_dataset(dataset: pd.DataFrame, disease_condition: list, perturbation_condition: list):

    new_df = pd.DataFrame(columns=dataset.columns)

    print(dataset.shape[0])

    for row in range(dataset.shape[0]):
        if dataset.iloc[row, 10] in disease_condition and dataset.iloc[row, 11] in perturbation_condition: #column 10 = disease,column 11 = perturbation
            current_row = dataset.loc[row]
            new_df = pd.concat([new_df, current_row.to_frame().T], ignore_index=True)


    total_cells = 0
    for row in range(new_df.shape[0]):
        print(f"Disease: {new_df.loc[row, "disease"]}. Perturbation: {new_df.loc[row, "perturbation"]}. Cells: {new_df.loc[row,"obs_count"]}")
        total_cells += new_df.loc[row, "obs_count"]

    print(f"The reduced dataframe contains {new_df.shape[0]} datasets with {total_cells} total cells")


    return new_df

# get a dataframe that has only datasets which are one of the relevant cancers and have no perturbation
potential_data = reduce_dataset(df, relevant_valid_disease_types, no_perturbation)

# 36 datasets for full filter (230K cells) 585 datasets for just cancer filter (4.3M cells)
# switch back to cancerSCEM + tabula sapiens
# RIP 12hours, at least I know scbasecount quite well now




