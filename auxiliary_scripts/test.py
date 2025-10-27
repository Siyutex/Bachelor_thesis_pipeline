import json
import scanpy as sc

adata = sc.read_h5ad()







"""
Example SLURM script

#!/bin/bash
#SBATCH --job-name=myjob
#SBATCH --output=myjob.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G

# Copy input files to local scratch
cp /home/username/myproject/data/input.csv $TMPDIR

# Run your Python script
python myscript.py --input $TMPDIR/input.csv --output $TMPDIR/result.csv

# Copy results back
cp $TMPDIR/result.csv /home/username/myproject/results/
"""





