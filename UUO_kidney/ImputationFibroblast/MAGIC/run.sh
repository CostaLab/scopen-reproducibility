#!//bin/bash


### Job name
#SBATCH -J run_magic
#SBATCH -e ./run_magic.txt
#SBATCH -o ./run_magic.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=180G -c 4


source ~/.bashrc
conda activate r-4.0.3

start=$(date +'%s')

time Rscript ./magic.R --input=../../Signac/data/Fibroblast/PeakMatrix.Rds --output=./MAGIC.txt

echo "It took $(($(date +'%s') - $start)) seconds"
