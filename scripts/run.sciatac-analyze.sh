#!/bin/bash

#
# Current date and time.
#
NOW=`date '+%Y%m%d.%H%M%S'`

#
# Path to the Nextflow processing run configuration file.
#
PWD=`pwd`
CONFIG_FILE="$PWD/experiment.config"

#
# Nextflow executable and pipeline script locations.
#
NEXTFLOW="/net/gs/vol1/home/bge/bin/nextflow"
NF_ANALYZE="/net/gs/vol1/home/bge/git/bbi-sciatac-analyze/main.nf"

#
# Get the path to the analyze output directory from
# the configuration file and set the Nextflow work
# directory to be in the analyze output directory.
# Notes:
#   o  bbi-sciatac-analyze/main.nf also defines ANALYZE_DIR as ${params.output_dir}/analyze_out so
#      changing the ANALYZE_DIR here requires a change in main.nf as well.
#   o  bbi-sciatac-analyze/main.nf defines DEMUX_DIR as ${params.output_dir}/demux_out so changes
#      to DEMUX_DIR requires a change in main.nf as well.
#
OUTPUT_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.output_dir"){print$2}}' | sed 's/"//g'`
ANALYZE_DIR="$OUTPUT_DIR/analyze_out"
WORK_DIR="$ANALYZE_DIR/work"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
REPORT_FIL="$ANALYZE_DIR/run_reports/analyze.report.html"
TRACE_FIL="$ANALYZE_DIR/run_reports/analyze.trace.tsv"
TIMELINE_FIL="$ANALYZE_DIR/run_reports/analyze.timeline.html"

#
# Nextflow run parameters.
#
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"

mkdir -p $ANALYZE_DIR/run_reports
pushd $ANALYZE_DIR

date > run_reports/run_start.${NOW}.txt

#
# Run Nextflow sci-ATAC analyze pipeline.
#
$NEXTFLOW run $NF_ANALYZE $PARS

date > run_reports/run_finish.${NOW}.txt

popd
