#!/bin/bash

### Job name
#SBATCH -J conformation
#SBATCH -e ./conformation.txt
#SBATCH -o ./conformation.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 3:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=100G -c 10

source ~/.bashrc
conda activate r-cicero

Rscript compute_aupr.R
