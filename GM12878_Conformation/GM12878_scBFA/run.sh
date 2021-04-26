#!/bin/bash

### Job name
#SBATCH -J conformation
#SBATCH -e ./conformation.txt
#SBATCH -o ./conformation.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 123:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=100G -c 10

source ~/.bashrc
conda activate r-4.0.3

Rscript /data/scATA/SingleCellOpenChromatin/local/ATAC/scripts/Imputation/scBFA.R --input ./GSM2970932_sciATAC_GM12878_matrix.txt --output ./GM12878.txt
#python /home/rs619065/scopen/scripts/convert2vector.py GM12878.txt GM12878_imputed.txt
#conda activate r-cicero
#Rscript RunCicero.R
