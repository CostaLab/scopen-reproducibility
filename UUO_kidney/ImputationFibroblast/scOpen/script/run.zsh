#!/bin/bash

source ~/.bashrc
conda activate py-3.6

start=$(date +'%s')

time scopen --input $1 --input_format 10X --output_dir $2 --output_prefix $3 --output_format dense --verbose 0 --nc 6 --estimate_rank 

echo "It took $(($(date +'%s') - $start)) seconds"
