- cell ids (obs_names) are consistent across RUNs (I tested this with compare_present_cell.py)
    -> I can compare what cells are labelled as transitional by obs_names

- there are only 61 cells which are consistently labelled as transitional acress RUNS 3.5, 3,6, 4 (also checked with compare_present_cells.py)
	- after makug pipeline deterministic, and subsampling 70% of cells for 5 runs of tree:
		12% of cells are transitional in all runs they appear in
		20% of cells are inconsistently labelled as transitional 
		(which is much better than the 1.4% vs 30% the where (wrongly) observed when considering
		 P(transitional) instead of P(transitioanl | sampled))

## cells with different expression across runs? ##
reduced adatas (X_scANVi_corrected_cnv) (3.5 3.6 4.0):
2467 cells with different expression in diff adatas 

aggregated with slightly diff mito cutoff:
0 cells with diff expr in diff adatas 
-> issue must happen after aggregation

reduced adatas (X_scANVI_corrected) (3.5 3.6 4.0):
2467 
-> same as cnv, so cnv algorithm must work deterministically -> must happen before cnv

-> must happen in batch correction (it was batch correction random sampling from ZINB + training new models every run)


## DOES STEP X PRODUCE CONSISTENT RESULTS? ## 
(UNDER THE CONDITION THAT INPUT AND HYPERPARAMETERS ARE CONSITENT!)

checks:
output across 3 runs needs to be like this:
	same cell ids present in all outputs
	same cell id needs to have same expression in all outputs for relevant layer
	same cell id needs to have same obs annotations in all outputs
	same gene ids present in all outputs
	same gene id needs to have same var annotations in all outputs

+ preprocessing.py: 
    -> passes all checks

+ Cell_type_annotation.py:
	-> passes all checks if same pretrained model used for all runs (same model used across all files)
		- suboptimal because different files might be from diff batches so no one sizes fits all
	-> fails assert_no_inconsistent_obs_annotations if a new model is trained per run (new model trained per file)
	-> passes all checks if same set of pretrained models is used for all runs (one pretrained model per file (based on consistent file names and contents))

- batch_aggregation.py:
	-> passes all checks

- Batch_correction.py:
	-> fails assert_no_inconsistent_expression, because get_normalized_expression draws a random sample (use a fixed seed for the RNG)
	-> passes all checks, for a set seed (69) with scvi.settings.seed = seed, expression matches for X_scANVI_corrected

- infer_CNV.py:
	-> passes all checks, expression matches for X_scANVI_corrected_cnv

- sc_malignant_finder.py
	-> passes all checks, expression matches for X_scANVI_corrected_cnv

- phylogenetic_tree.py
	-> only fails assert_no_inconsistent_obs_annotations, most likely because tree building is non
	deterministic
	-> passes all checks, for a set seed for the numpy rng for adding jitter to distance metric


- matrix_isolation_HVGs.py
	-> passes all checks, expression matches for log1p
	