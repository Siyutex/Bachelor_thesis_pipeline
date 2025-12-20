# This script takes scRNAseq data (shouled be preprocessed, aggregatedand batch corrected) and
# inferes the edges of a GRN from it using grnboost2
# the number of edges can be limited by max_n_edges or min_importance_score

# ye olde bootstrappin' plan:
"""
- Subsample N cells, with replication for an N cell expression matrix (duplicates can occur, p of a cell occuring = 1-(1-1/N)^N ~ 0.63)
- run grnboost2 with those cells
- record top 10 edges per target in a sparse matrix (BxE, bootstraps times edges) with 1 as the entry (1 = occurred in this bootstrap run)
- repeat continuously
- for each edge, compute stdev across bootstraps, then average across all edges
- track delta_total_stdev, when it gets below 5% of total_stdev, stop bootstrapping
- only take the edges with low stdev and mean close to 1 (= on average, they are more likely to occur)
- export edges as list of tuples [(target <- regulator), ...]
"""

import scanpy as sc
import pandas as pd
import numpy as np
import helper_functions as hf
import os
import dask
dask.config.set({"dataframe.query-planning": False}) # need to not get import errors
from arboreto.algo import grnboost2
import json


class StabilityTracker:
    def __init__(self, convergence_threshold=0.05, top_n=5, min_runs=20, min_stability=0.5):
        self.convergence_threshold = convergence_threshold # convergence limit; stops when average edge stability (stdev) changes by less than tol from one run to next
        self.top_n = top_n # top n most stable / frequent regulators per target
        self.min_runs = min_runs # minimum amount of runs before converging (it could otherwise happend that stdevs randomly fall close to eachother and premture false convergence occurs)
        self.min_stability = min_stability  # putative regulators need to occur at least in this percentage of rens eg 0.5 = 50% of runs
        self.run_count = 0
        self.converged = False
        
        # Track counts of (Target, TF) pairs. Index: MultiIndex(target, TF) (need to initialize wih multiindex labels)
        self.edge_counts = pd.Series(
            dtype=float, 
            index=pd.MultiIndex.from_tuples([], names=["target", "TF"])
        )
        self.prev_overall_stdev = None

    def add_run(self, series):
        self.run_count += 1
        
        # 1. Convert current series to a flat set of edges
        records = []
        for target, tf_list in series.items():
            for tf in tf_list:
                records.append((target, tf)) # list of tuples of (target, TF) -> basically what edges were relevant in this run

        # convert to pandas multiindex
        current_edges = pd.MultiIndex.from_tuples(records, names=["target", "TF"])  

        # 2. Update edge_counts (Running Sum)
        # Using .add with fill_value=0 keeps the index growing as new edges appear
        new_hits = pd.Series(1, index=current_edges) # we use the current edges multiindex to index the series with values 1 (the edge existed in the run)
        self.edge_counts = self.edge_counts.add(new_hits, fill_value=0) # new_hits added to running sum for eac edge (if the edge existed, you get +1 vote for this edge, if it didn't +0 (fill value))

        # 3. Calculate Stats on the fly
        # Mean = count / total_runs
        means = self.edge_counts / self.run_count # pandas series, index = edges (multiindex), value = current relative frequency of the edge across runs
        
        # Binary Stdev formula: sqrt(p * (1-p))
        # This is much faster than calculating std() on a huge 0/1 matrix
        stdevs = np.sqrt(means * (1 - means))
        overall_stdev = stdevs.mean()

        # 4. Convergence Check
        if self.prev_overall_stdev is not None:
            rel_change = abs((overall_stdev - self.prev_overall_stdev) / self.prev_overall_stdev)
            # Only converge after min_runs to avoid early noise
            if rel_change < self.convergence_threshold and self.run_count >= self.min_runs:
                self.converged = True
        else:
            rel_change = np.inf # in the first run, where there is no prev_overall_stdev, there is no sensical relative change

        self.prev_overall_stdev = overall_stdev
        
        vprint(f"Run {self.run_count} | Stdev: {overall_stdev:.4f} | Change: {rel_change:.4f}")
        return self.converged

    def get_final_selection(self):
        if self.edge_counts.empty:
            raise ValueError("No runs added yet.")

        # 1. Calculate final frequencies (0.0 to 1.0)
        means = self.edge_counts / self.run_count
        
        # 2. Filter by minimum stability FIRST
        # This removes any TF that appeared in fewer than min_stability % of runs
        stable_means = means[means >= self.min_stability]

        # 3. Select top N from the remaining stable pool

        # grouby(level = 0) splits the series into buckets based in the first entry in the multiindex (so target in our case)
        # then we only keep the top_n TFs in each bucket based on the relative frequency
        # result is still a pandas series but with only top_n TFs per target
        final_selection = ( 
            stable_means.groupby(level="target", group_keys=False)
                        .nlargest(self.top_n)
        )

        # 4. Convert to dictionary
        result = {}
        # grouby yields tuples. the first entry is the thing that should be grouped by (here the target), the second one a sub series of the original series that has that target in the multiindex
        # we then access this sub series' index and get the values at the "TF" level, turn to list
        # this results in a dictionary with key = target and value = list of TFs
        for target, group_series in final_selection.groupby(level=0): 
            result[target] = group_series.index.get_level_values("TF").tolist()
            
        return result
    

def get_tf_ensg_list(adata, tf_file_path):
    """
    Maps a list of TF gene symbols to ENSG IDs using the mapping in adata.var.
    
    Parameters:
    -----------
    adata : sc.AnnData
        The AnnData object where index is ENSG and var['gene_symbols'] is Symbols.
    tf_file_path : str
        Path to the downloaded TF text file.
    """
    # 1. Load the TF symbols from the file
    with open(tf_file_path, 'r') as f:
        # Strip whitespace and ignore empty lines or headers
        external_tfs = {line.strip() for line in f if line.strip()}

    # 2. Create a mapping dictionary: {Symbol: ENSG}
    # We swap the index and the column to make lookups fast
    symbol_to_ensg = pd.Series(adata.var_names, index=adata.var['gene_symbols']).to_dict()

    # 3. Intersection and Mapping
    # We only keep TFs that are actually present in your dataset's var_names
    tf_ensg_list = [
        symbol_to_ensg[sym] 
        for sym in external_tfs 
        if sym in symbol_to_ensg
    ]

    vprint(f"Loaded {len(external_tfs)} TFs from file.")
    vprint(f"Mapped {len(tf_ensg_list)} TFs to ENSG IDs present in adata.")
    
    return tf_ensg_list





if __name__ == "__main__":
   
    # get cmd args
    input_data_file, output_data_dir, layer, n_nodes, tf_list_file, convergence_threshold, top_n_regulators, min_runs, min_stability, input_prefix, output_prefix, verbose = hf.import_cmd_args(12)
    vprint = hf.make_vprint(verbose)

    # intialize RNG
    np.random.seed(42)
    rng = np.random.default_rng(seed=42)

    # import data from h5ad file
    vprint(f"Importing data from {input_data_file}")
    if layer == "X":
        adata = sc.read_h5ad(input_data_file)
    else:
        adata = hf.matrix_to_anndata(sc.read_h5ad(input_data_file), layer)

    # check if obs and var names are unique
    if len(adata.obs_names.unique()) != adata.n_obs:
        raise ValueError("Obs names are not unique, obs names should be made unique in aggregate_batches during batch correction.")
    if len(adata.var_names.unique()) != adata.n_vars:
        raise ValueError("Var names are not unique, var names should be unique be default (reference genome should not contain duplicates).")

    # if n_nodes is not None, subsample genes randomly (used for debugging, not recommended for real data)
    if n_nodes is not None:
        vprint(f"Limiting GRN scope to {n_nodes} nodes...")
        idx = rng.choice(adata.X.shape[0], size=n_nodes, replace=False, shuffle=False)
        adata = adata[:, idx]

    # get TF list
    tf_ensg_list = get_tf_ensg_list(adata, tf_list_file) if tf_list_file is not None else "all"

    tracker = StabilityTracker(convergence_threshold=convergence_threshold, top_n=top_n_regulators, min_runs=min_runs, min_stability=min_stability)
    while tracker.converged == False:
        # subsrample the dataframe randomly (adata.shape[0] random samples WITH replacement = bootstrapping)
        idx = rng.choice(adata.X.shape[0], size=adata.X.shape[0], replace=True, shuffle=False)

        # use those indices to pick rows
        subsample = adata[idx, :]
        vprint(f"Adata shape after subsampling: {subsample.X.shape}")

        # assign dataframe with var names as column names
        subsample_df = pd.DataFrame(subsample.X.toarray(), index=subsample.obs_names, columns=subsample.var_names)
        vprint("Dataframe shape: ", subsample_df.shape)

        # run GRN inference (^2 compute time, 2324 genes take 2:40 minutes, cells do not seem to affect runtime)
        vprint("Running GRN inference...")
        grn = grnboost2(subsample_df, verbose=verbose, tf_names=tf_ensg_list, seed=42)

        # Calculate the 95th percentile threshold of the 'importance' column
        importance_threshold = grn['importance'].quantile(0.95)
        vprint(f"95th percentile importance threshold: {importance_threshold:.4f}")

        # Filter for high-confidence edges
        # We also sort by importance to ensure the lists are ordered by strength
        top_grn = (grn[grn['importance'] >= importance_threshold]
                   .sort_values("importance", ascending=False))

        # Create the series of TF lists for the tracker
        # groupby target -> splits DF into 1 bucket per unique itme in target
        # ["TF"] -> makes sure that only the TF column gets to be values in that bucket
        # apply(list) -> turns each bucket into a list
        # TF_list = pandas series with index = target gene and values = >95th importance score percentile regulators
        TF_list = top_grn.groupby("target")["TF"].apply(list) 

        # update tracker with result from current bootsrapping run
        tracker.add_run(TF_list)

    
    if tracker.converged == True:
        final_selection = tracker.get_final_selection()
        output_file_path = os.path.join(output_data_dir, f"{output_prefix}_{os.path.basename(input_data_file).removeprefix(input_prefix).removesuffix('.h5ad')}.json")
    
        with open(output_file_path, "w") as f:
            json.dump(final_selection, f)

    # send output to executor
    print(f"Output: {output_file_path}")