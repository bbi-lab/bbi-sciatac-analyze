#!/bin/bash


# sample_name=$1
# matrix_path=$2
# umi_cutoff=$3

run_name='ATAC2-001'
sample_name='ATAC.BAT1'

run_name='ATAC3-002'
sample_name='E14.5'


matrix_path="../$run_name/make_matrices/${sample_name}-peak_matrix.mtx.gz"
umi_cutoff=60

Rscript run_scrublet.R --sample_name=$sample_name --matrix_path=$matrix_path --umi_cutoff=$umi_cutoff
run_scrublet.py --mat="${sample_name}_to_scrublet.mtx" --key=$sample_name
