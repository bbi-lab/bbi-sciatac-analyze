#!/bin/bash


function make_aligner_index()
{
  echo "Build bowtie2 index files ..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  bowtie2-build $FASTA_FILTERED $INDEX_PREFIX 2>&1 | tee -a ${LOG_FILE}

  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}
