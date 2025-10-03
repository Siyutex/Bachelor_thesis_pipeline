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
import sys
import helper_functions as hf
import os
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
import json

# get cmd args
input_data_file, output_data_dir, verbose = hf.import_cmd_args(3)
vprint = hf.make_vprint(verbose)


class StabilityTracker:
    """
    Class for tracking stability of TFs inferred by grnboost2 + utility functions.
    """
    def __init__(self, tol=0.05, top_n=10):
        """
        tol   = relative change threshold for convergence
        top_n = number of most stable TFs to export per target
        """
        self.tol = tol
        self.top_n = top_n
        self.stability_df = pd.DataFrame(columns=[], index=pd.MultiIndex.from_tuples([], names=["target", "TF"]))  # MultiIndex (target, TF), cols = bootstrap runs
        self.prev_overall_stdev = None
        self.run_count = 0
        self.means = None
        self.stdevs = None
        self.overall_stdev = None
        self.converged = False

    def add_run(self, series):
        """
        Add one bootstrap run.
        series: pandas Series, index=targets, values=list of TFs (top TFs targetting that target).
        """
        # Build run dataframe
        records = []
        for target, tf_list in series.items():
            for tf in tf_list:
                records.append((target, tf, 1))
        run_df = pd.DataFrame(records, columns=["target", "TF", self.run_count])
        vprint(f"Dataframe for run {self.run_count} :\n{run_df.head(5)}")
        run_df = run_df.set_index(["target", "TF"])

        # Merge into stability_df
        self.stability_df = self.stability_df.join(run_df, how="outer").fillna(0)

        # Update stats
        self.means = self.stability_df.mean(axis=1)
        self.stdevs = self.stability_df.std(axis=1)
        self.overall_stdev = self.stdevs.mean()

        # print run stats
        vprint(f"Run {tracker.run_count}, overall stdev: {tracker.overall_stdev}")
        vprint(f"Previous overall stdev: {tracker.prev_overall_stdev}")

        # Check convergence
        rel_change = None
        if self.prev_overall_stdev is not None and self.overall_stdev is not None:
            rel_change = abs((self.overall_stdev - self.prev_overall_stdev) / self.prev_overall_stdev)
            if rel_change < self.tol:
                self.converged = True
        vprint(f"Relative change in stdev: {rel_change}")
        vprint(f"Converged: {tracker.converged}")

        self.prev_overall_stdev = self.overall_stdev
        self.run_count += 1

        return self.converged

    def get_final_selection(self):
        """
        After convergence, pick top_n most stable TFs per target (highest mean).
        Stdev is a bad measure because it can also select consistently absent edges.
        Higher mean is directly proportional to lower stdev, sow mean as the proxy
        for stability should be sufficient.
        Returns dict[target] = list of TFs
        """
        if self.means is None:
            raise ValueError("No runs added yet.")

        final_selection = (
            self.means.groupby("target")
                       .nlargest(self.top_n)
                       .reset_index(level=0, drop=True)
        )

        result = {}
        for target, group in final_selection.groupby("target"):
            result[target] = group.index.get_level_values("TF").tolist()
        return result

    def export_json(self, filepath="stable_tfs.json"):
        result = self.get_final_selection()
        with open(filepath, "w") as f:
            json.dump(result, f, indent=2)
        return filepath



if __name__ == "__main__":
    # consts
    PREFIX = "batch_corrected_HVG_"
    SUFFIX = ".h5ad"




    # import data from batch corrected h5ad file
    vprint(f"Importing data from {input_data_file}")
    adata = sc.read_h5ad(input_data_file)

    # check if obs and var names are unique
    if len(adata.obs_names.unique()) != adata.n_obs:
        raise ValueError("Obs names are not unique, obs names should be made unique in aggregate_batches during batch correction.")
    if len(adata.var_names.unique()) != adata.n_vars:
        raise ValueError("Var names are not unique, var names should be unique be default (reference genome should not contain duplicates).")

    # isolate ductal cells
    vprint("Isolating ductal cells...")
    ductal_cells = adata.obs["cell_type"] == "ductal_cell"
    adata = adata[ductal_cells, :]

    tracker = StabilityTracker(tol=0.05, top_n=5)
    while tracker.converged == False:
        # subsrample the dataframe randomly
        rng = np.random.default_rng()
        idx = rng.choice(adata.X.shape[0], size=adata.X.shape[0], replace=True, shuffle=False)
        # print(idx) # debugging
        # print(len(list(set(idx)))) # debugging
        idx = list(set(idx))

        # use those indices to pick rows
        subsample = adata[idx, :]
        print(f"Adata shape after subsampling: {subsample.X.shape}")

        # subsample = subsample[:, :100] # debugging

        # assign dataframe with var names as column names
        subsample_df = pd.DataFrame(subsample.X.toarray(), index=subsample.obs_names, columns=subsample.var_names)
        print("Dataframe shape: ", subsample_df.shape)

        # assign client
        # client = Client(LocalCluster())

        # run GRN inference (^2 compute time, 2324 genes take 2:40 minutes, cells do not seem to affect runtime)
        vprint("Running GRN inference...")
        grn = grnboost2(subsample_df, verbose=verbose, tf_names="all")

        # get a pandas series with top 10 edges per target for current run
        # print(type(grn))
        # print(grn.groupby("target").head(10))
        TF_list = grn.groupby("target")["TF"].apply(list) # groupby returns one dataframe per unique entry in the column passed as an argument, apply then runs the passed functions on each of those dataframes, and returns a series (index = dataframe name, value = returned value of passed function (in this case, a list of all the entries in that dataframe))
        # print(type(TF_list))
        # print(TF_list)

        # update tracker with result from current bootsrapping run
        tracker.add_run(TF_list)


    
    if tracker.converged == True:
        # export json with dict of final TF selection
        output_file_path = tracker.export_json(filepath=os.path.join(output_data_dir, f"Stable_edges_{os.path.basename(input_data_file).removeprefix(PREFIX).removesuffix(SUFFIX)}.json"))        

    # send output to executor
    print(f"Output: {output_file_path}")

    sys.exit(0)