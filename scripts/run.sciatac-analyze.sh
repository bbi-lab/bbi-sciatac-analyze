#!/bin/bash

#
# Edit this file to set the NEXTFLOW variable to the location
# of the Nextflow executable, and NF_MAIN to the location of
# the bbi-sciatac-demux Nextflow pipeline script, main.nf.
#

# 
# Copy this file to the directory where you will run the pipeline.
# The experiment.config file must be in the directory where you
# run this script.
#

#
# Nextflow executable and pipeline script locations.
#
NEXTFLOW="<path_to_nextflow_program>"
NF_MAIN="<path_to_bbi-sciatac-analyze_repository/main.nf"

#
# Current date and time.
#
NOW=`date '+%Y%m%d_%H%M%S'`

#
# Path to the Nextflow processing run configuration file.
#
PWD=`pwd`
CONFIG_FILE="$PWD/experiment.config"

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
# REPORT_FIL="$ANALYZE_DIR/run_reports/analyze.report.${NOW}.html"
TRACE_FIL="$ANALYZE_DIR/run_reports/analyze.trace.${NOW}.tsv"
# TIMELINE_FIL="$ANALYZE_DIR/run_reports/analyze.timeline.${NOW}.html"

#
# Nextflow run parameters.
#
# PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FIL -resume"

mkdir -p $ANALYZE_DIR/run_reports
pushd $ANALYZE_DIR

date > run_reports/run_start.${NOW}.txt

#
# Run Nextflow sci-ATAC analyze pipeline.
#
$NEXTFLOW run $NF_MAIN $PARS

date > run_reports/run_finish.${NOW}.txt

popd
