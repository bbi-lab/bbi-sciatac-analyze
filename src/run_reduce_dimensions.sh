#!/bin/bash

lsample='BAT100'
root_dir='/net/bbi/vol1/data/sciATACseq/nextseq_runs/ATAC3-002-a/analyze_out'
script_dir='/net/gs/vol1/home/bge/git/bbi-sciatac-analyze/src'

function run_reduce_dimensions()
{
  local l_sample_name=${1}
  local l_root_dir=${2}

  local l_mat_dir="${l_root_dir}/${l_sample_name}/make_matrices"
  local l_count_dir="${l_root_dir}/${l_sample_name}/count_report"
  local l_umap_plot="${l_root_dir}/${l_sample_name}/reduce_dimension/${l_sample_name}-umap_plot.png"
  local l_cds_file="${l_root_dir}/${l_sample_name}/reduce_dimension/${l_sample_name}-monocle3_cds.rds"

echo "mat_dir: $l_mat_dir"
echo "count_dir: $l_count_dir"
echo "umap_plot: $l_umap_plot"
echo "cds_file: $l_cds_file"


#  black_list_file='mm10-blacklist.v2.sorted.bed'
  
  
  umi_cutoff=100
  frip_cutoff=0.1
  frit_cutoff=0.05
  num_lsi_dimensions=75
  
  
  $script_dir/run_scrublet.py --sample_name=$l_sample_name --mat=${l_mat_dir}/${l_sample_name}-peak_matrix.mtx.gz --umi_cutoff=$umi_cutoff
  
  Rscript $script_dir/reduce_dimensions.R --mat_dir=$l_mat_dir \
                                          --count_dir=$l_count_dir \
                                          --sample_name=$l_sample_name \
                                          --umi_cutoff=$umi_cutoff \
                                          --frip_cutoff=$frip_cutoff \
                                          --frit_cutoff=$frit_cutoff \
                                          --doublet_predict \
                                          --doublet_predict_top_ntile=0.1 \
                                          --num_lsi_dimensions=$num_lsi_dimensions \
                                          --umap_plot=$l_umap_plot \
                                          --cds_file=$l_cds_file

#                                          --black_list_file=$black_list_file \
#                                          --lsi_coords_file=$lsi_coords_file \
#                                          --umap_coords_file=$umap_coords_file \
}


for sample in $lsample
do
  run_reduce_dimensions $sample $root_dir
done

