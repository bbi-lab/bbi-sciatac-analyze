#!/bin/bash

#
# Path to the Nextflow processing run configuration file.
#
CONFIG_FILE="/xxx/params.config"

#
# Nextflow executable and pipeline script locations.
#
NEXTFLOW="/xxx/bin/nextflow"
NF_ANALYZE="/xxx/bbi-sciatac-analyze/main.nf"

#
# Get the path to the analyze output directory from
# the configuration file and set the Nextflow work
# directory to be in the analyze output directory.
#
ANALYZE_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.analyze_dir"){print$2}}' | sed 's/"//g'`
WORK_DIR="$ANALYZE_DIR/work"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
REPORT_FIL="$ANALYZE_DIR/analyze.report.html"
TRACE_FIL="$ANALYZE_DIR/analyze.trace.tsv"
TIMELINE_FIL="$ANALYZE_DIR/analyze.timeline.html"

#
# Nextflow run parameters.
# Note: I include -resume for convenience. I does not affect the initial run.
#
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"

mkdir -p $ANALYZE_DIR
pushd $ANALYZE_DIR

date > ./run_start.txt

#
# Run Nextflow sci-ATAC analyze pipeline.
#
$NEXTFLOW run $NF_ANALYZE $PARS

date > run_finish.txt

popd
