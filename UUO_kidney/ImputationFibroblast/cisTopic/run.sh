#!/bin/bash

### Job name
#SBATCH -J run_cisTopic
#SBATCH -e ./run_cisTopic.txt
#SBATCH -o ./run_cisTopic.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=180G -c 4

source ~/.bashrc
conda activate r-4.0.3

start=$(date +'%s')

time Rscript ./cisTopic.R --input=../../Signac/data/Fibroblast/PeakMatrix.Rds --output=./

echo "It took $(($(date +'%s') - $start)) seconds"
