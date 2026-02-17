# Bachelor\_thesis\_pipeline

This repository contains the code associated with the thesis "Implementation and Validation of a Computational Pipeline for Reversion Switch Identification in Primary PDAC". A guide to initializing the necessary miniconda3 environment can be found under "pipeline_environment". Model weights can be found under models_PDAC (for the PDAC dataset referenced in the thesis) and models_shin (for the original REVERT dataset referenced in the thesis). All scripts but the pipeline_executor contain the "default parameters" and default seeds referenced in the thesis. The pipeline executor contains commented out blocks of code in the code execution section at the end that contain the default parameters for the subprocess scripts. A list of libraries directly relevant for imports can be found in the file "used_libraries.txt"

IMPORTANT INFO: 
-there is a bug in arboreto 0.1.6 that yields the error "must supply at least on delayed object"; to fix go into arboreto.core.create\_graph and comment out the line that says "all\_meta\_df = from\_delayed(delayed\_meta\_dfs, meta=\_META\_SCHEMA)"

