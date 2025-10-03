# takes edge set from edge rule inference 
# (pseudotime binning, several runs)
# use scikit learn to get a random forest for each node + regulators, choose best decision tree
# turn decision tree into boolean function
# espresso minimization
# export rule set as csv, then use for pyboolnet simulation / analysis

import sklearn
import pyeda
import pandas as pd
import json
import scanpy as sc
import helper_functions as hf

# import cmd args
expression_matrix_file, output_data_dir, verbose, expression_matrix_file = hf.import_cmd_args(4)
vprint = hf.make_vprint(verbose)

# load edge set from json (list of tuples)
with open("edge_set.json", "r") as f:
    edge_set = json.load(f)
edge_set = edge_set["edge_set"]

# load expression matrix
adata = sc.read_h5ad(expression_matrix_file)







