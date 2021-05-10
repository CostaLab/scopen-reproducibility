#!/usr/local_rwth/bin/zsh

### Job name
#SBATCH -J run_SCALE
#SBATCH -e ./run_SCALE.txt
#SBATCH -o ./run_SCALE.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 1:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=180G -c 4 --partition c18g --gres gpu:1


source ~/.zshrc
conda activate r-4.0.3


start=$(date +'%s')

time python ~/SCALE/SCALE.py -d ../../Signac/data/Fibroblast/PeakMatrix.txt -o ./ --impute --seed 42 --min_peaks 0 -x 0 \
--high 1.0

echo "It took $(($(date +'%s') - $start)) seconds"
