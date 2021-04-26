#!/usr/bin/env zsh

### Job name
#SBATCH -J conformation
#SBATCH -e ./conformation.txt
#SBATCH -o ./conformation.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 3:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=60G -c 10 --partition c18g --gres gpu:1

module load cuda/100
module load cudnn/7.4

source ~/.zshrc
conda activate r-4.0.3

python ~/SCALE/SCALE.py -d GSM2970932_sciATAC_GM12878_matrix.txt -o ./ --impute --seed 42 --lr 0.0002 --min_peaks 0 -x 0

