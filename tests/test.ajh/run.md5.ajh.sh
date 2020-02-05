#!/bin/bash

dst="$1"

lfil='md5.ajh.align_reads.sh md5.ajh.merge_bams.sh md5.ajh.call_cells.sh md5.ajh.call_peaks.sh md5.ajh.count_report.sh md5.ajh.get_unique_fragments.sh md5.ajh.tss_region_coverage.sh md5.ajh.make_matrices.sh md5.ajh.summarize_cell_calls.sh'

function run_check()
{
  local nmscript="$1"
  test_dir=`echo $nmscript | sed 's/md5.ajh.//' | sed 's/.sh//'`
  if [ "$test_dir" == "tss_region_coverage" ]
  then
    test_dir="per_base_tss_region_coverage"
  fi

  pushd ${dst}/${test_dir}
  ${nmscript}
  popd
}

for fil in $lfil
do
  run_check $fil
done

