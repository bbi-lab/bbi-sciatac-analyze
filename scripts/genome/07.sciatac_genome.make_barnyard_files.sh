#!/bin/bash

function barnyard_edit_fastas()
{
  HUMAN_FASTA=`echo $HUMAN_FASTA_GZ | sed 's/\.gz//'`
  HUMAN_FASTA_FILTERED=${HUMAN_FASTA}.filtered
  HUMAN_FASTA_EDITED=${HUMAN_FASTA}.filtered.edited
  echo "Edit fasta file ${HUMAN_FASTA_FILTERED}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  sed 's/^>/>HUMAN_/' $HUMAN_FASTA_FILTERED > $HUMAN_FASTA_EDITED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Calculate fasta file sequence md5 checksums for ${HUMAN_FASTA_EDITED}..." | tee -a ${LOG_FILE}
  $MD5_SEQ $HUMAN_FASTA_EDITED > ${HUMAN_FASTA_EDITED}.md5_seq
  cat $HUMAN_FASTA_EDITED.md5_seq | tee -a ${LOG_FILE}

  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}


  MOUSE_FASTA=`echo $MOUSE_FASTA_GZ | sed 's/\.gz//'`
  MOUSE_FASTA_FILTERED=${MOUSE_FASTA}.filtered
  MOUSE_FASTA_EDITED=${MOUSE_FASTA}.filtered.edited
  echo "Edit fasta file ${MOUSE_FASTA_FILTERED}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  sed 's/^>/>MOUSE_/' $MOUSE_FASTA_FILTERED > $MOUSE_FASTA_EDITED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Calculate fasta file sequence md5 checksums for ${MOUSE_FASTA_EDITED}..." | tee -a ${LOG_FILE}
  $MD5_SEQ $MOUSE_FASTA_EDITED > ${MOUSE_FASTA_EDITED}.md5_seq
  cat $MOUSE_FASTA_EDITED.md5_seq | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function barnyard_edit_chromosome_sizes()
{
  HUMAN_CHROMOSOME_SIZES_FILE_EDITED=${HUMAN_CHROMOSOME_SIZES_FILE}.edited
  echo "Edit chromosome sizes file ${HUMAN_CHROMOSOME_SIZES_FILE}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  sed 's/^/HUMAN_/' $HUMAN_CHROMOSOME_SIZES_FILE > $HUMAN_CHROMOSOME_SIZES_FILE_EDITED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $HUMAN_CHROMOSOME_SIZES_FILE_EDITED | tee -a ${LOG_FILE}

  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}


  MOUSE_CHROMOSOME_SIZES_FILE_EDITED=${MOUSE_CHROMOSOME_SIZES_FILE}.edited
  echo "Edit chromosome sizes file ${MOUSE_CHROMOSOME_SIZES_FILE}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  sed 's/^/MOUSE_/' $MOUSE_CHROMOSOME_SIZES_FILE > $MOUSE_CHROMOSOME_SIZES_FILE_EDITED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $MOUSE_CHROMOSOME_SIZES_FILE_EDITED | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function barnyard_edit_whitelist_regions()
{
  HUMAN_WHITELIST_REGIONS_BED_EDITED=${HUMAN_WHITELIST_REGIONS_BED}.edited
  echo "Edit whitelist regions file ${HUMAN_WHITELIST_REGIONS_BED}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE} 
  sed 's/^/HUMAN_/' $HUMAN_WHITELIST_REGIONS_BED > $HUMAN_WHITELIST_REGIONS_BED_EDITED
  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $HUMAN_WHITELIST_REGIONS_BED_EDITED | tee -a ${LOG_FILE}


  MOUSE_WHITELIST_REGIONS_BED_EDITED=${MOUSE_WHITELIST_REGIONS_BED}.edited
  echo "Edit whitelist regions file ${MOUSE_WHITELIST_REGIONS_BED}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  sed 's/^/MOUSE_/' $MOUSE_WHITELIST_REGIONS_BED > $MOUSE_WHITELIST_REGIONS_BED_EDITED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $MOUSE_WHITELIST_REGIONS_BED_EDITED | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function barnyard_edit_tss_bed()
{
  HUMAN_TSS_BED_EDITED=${HUMAN_TSS_BED}.edited
  echo "Edit tss bed file ${HUMAN_TSS_BED}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE} 
  zcat ${HUMAN_TSS_BED}.gz | sed 's/^/HUMAN_/' > ${HUMAN_TSS_BED_EDITED}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count entries by sequence name in ${HUMAN_TSS_BED_EDITED}..." | tee -a ${LOG_FILE}
  cat ${HUMAN_TSS_BED_EDITED} | awk '{print$1}' | sort | uniq -c | sort -k 2,2V | awk '{printf( "%s\t%s\n", $1, $2 );}' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  MOUSE_TSS_BED_EDITED=${MOUSE_TSS_BED}.edited
  echo "Edit tss bed file ${MOUSE_TSS_BED}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  zcat ${MOUSE_TSS_BED}.gz | sed 's/^/MOUSE_/' > ${MOUSE_TSS_BED_EDITED}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count entries by sequence name in ${MOUSE_TSS_BED_EDITED}..." | tee -a ${LOG_FILE}
  cat ${MOUSE_TSS_BED_EDITED} | awk '{print$1}' | sort | uniq -c | sort -k 2,2V | awk '{printf( "%s\t%s\n", $1, $2 );}' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function barnyard_combine_files()
{
  BARNYARD_FASTA="barnyard.fa"
  echo "Combine human and mouse fasta files into ${BARNYARD_FASTA}..."
  cat $HUMAN_FASTA_EDITED $MOUSE_FASTA_EDITED > $BARNYARD_FASTA

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Calculate fasta file sequence md5 checksums for $BARNYARD_FASTA..." | tee -a ${LOG_FILE}
  $MD5_SEQ $BARNYARD_FASTA > ${BARNYARD_FASTA}.md5_seq
  cat ${BARNYARD_FASTA}.md5_seq | tee -a ${LOG_FILE}

  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  BARNYARD_CHROMOSOME_SIZES_FILE="barnyard.chromosome_sizes.txt"
  echo "Combine human and mouse chromosome sizes files into ${BARNYARD_CHROMOSOME_SIZES_FILE}..."
  cat  $HUMAN_CHROMOSOME_SIZES_FILE_EDITED $MOUSE_CHROMOSOME_SIZES_FILE_EDITED > $BARNYARD_CHROMOSOME_SIZES_FILE

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $BARNYARD_CHROMOSOME_SIZES_FILE | tee -a ${LOG_FILE}

  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  BARNYARD_WHITELIST_REGIONS_BED="barnyard.whitelist_regions.bed"
  echo "Combine human and mouse chromosome sizes files into ${BARNYARD_WHITELIST_REGIONS_BED}..."
  cat $HUMAN_WHITELIST_REGIONS_BED_EDITED $MOUSE_WHITELIST_REGIONS_BED_EDITED > $BARNYARD_WHITELIST_REGIONS_BED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $BARNYARD_WHITELIST_REGIONS_BED | tee -a ${LOG_FILE}

  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  BARNYARD_TSS_BED="barnyard.tss.bed"
  echo "Combine human and mouse tss files into ${BARNYARD_TSS_BED}..."
  cat ${HUMAN_TSS_BED_EDITED} ${MOUSE_TSS_BED_EDITED} | gzip > ${BARNYARD_TSS_BED}.gz

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count entries by sequence name in ${BARNYARD_TSS_BED}..." | tee -a ${LOG_FILE}
  zcat ${BARNYARD_TSS_BED}.gz | awk '{print$1}' | sort | uniq -c | sort -k 2,2V | awk '{printf( "%s\t%s\n", $1, $2 );}' | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo "Done." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function barnyard_make_aligner_index()
{
  echo "Build bowtie2 index files for ${BARNYARD_FASTA}..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  bowtie2-build $BARNYARD_FASTA $INDEX_PREFIX 2>&1 | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function barnyard_remove_unnecessary_files()
{
  echo "hi"
}
