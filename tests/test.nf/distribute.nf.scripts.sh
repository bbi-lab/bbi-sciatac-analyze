#!/bin/bash

dst="$1"

cp md5.nf.align_reads.sh ${dst}/align_reads
cp md5.nf.call_cells.sh ${dst}/call_cells
cp md5.nf.call_peaks.sh ${dst}/call_peaks
cp md5.nf.count_report.sh ${dst}/count_report
cp md5.nf.get_unique_fragments.sh ${dst}/get_unique_fragments
cp md5.nf.make_matrices.sh ${dst}/make_matrices
cp md5.nf.merge_bams.sh ${dst}/merge_bams
cp md5.nf.summarize_cell_calls.sh ${dst}/summarize_cell_calls
cp md5.nf.tss_region_coverage.sh ${dst}/per_base_tss_region_coverage

