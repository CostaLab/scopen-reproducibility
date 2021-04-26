#!/usr/bin/env zsh

### Job name
#SBATCH -J conformation
#SBATCH -e ./conformation.txt
#SBATCH -o ./conformation.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 3:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=100G -A rwth0233 -c 10

source ~/.zshrc
conda activate r-4.0.3

#scopen --input GSM2970932_sciATAC_GM12878_matrix.txt --output_dir ./ --output_prefix GM12878 --estimate_rank --nc 10
#python /home/rs619065/SingleCellOpenChromatin/utils/convert2vector.py GM12878.txt GM12878_imputed.txt
Rscript RunCicero.R
