#!/bin/bash


#
# Choose names of Ensembl REF-type sequences.
#
function sequences_to_keep_ref()
{
  echo "Keep the REF sequences for read alignments..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
  cat $FASTA.headers | awk '{if($4=="REF"){print$1}}' | sed 's/^>//' > $FINAL_IDS_FILE

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $FINAL_IDS_FILE | tee -a ${LOG_FILE}

  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


#
# Choose sequences named in the environment variable SEQUENCES_TO_KEEP_ALIGNER...
#
function sequences_to_keep_named()
{
  echo "Keep the named sequences for read alignments" | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
  echo "$SEQUENCES_TO_KEEP_ALIGNER" | sed 's/[ ][ ]*/\n/g' > $FINAL_IDS_FILE

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $FINAL_IDS_FILE | tee -a ${LOG_FILE}

  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


#
# Extract the specified sequences into a new fasta file.
#
function filter_fasta_file()
{ 
  echo "Extract selected sequences from the FASTA file..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  /net/bbi/vol1/data/src/sequtil/fasta_getseqs $FINAL_IDS_FILE $FASTA > $FASTA_FILTERED
  $FASTA_GETSEQS $FINAL_IDS_FILE $FASTA > $FASTA_FILTERED

  echo "Calculate filtered FASTA file sequence md5 checksums..." | tee -a ${LOG_FILE}
  $MD5_SEQ $FASTA_FILTERED > $FASTA_FILTERED.md5_seq

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Compare the md5 checksums for the downloaded and filtered FASTA files..." | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
  sort -k1,1 $FASTA.md5_seq > $FASTA.md5_seq.sort
  sort -k1,1 $FASTA_FILTERED.md5_seq > $FASTA_FILTERED.md5_seq.sort
  join -1 1 -2 1 $FASTA.md5_seq.sort $FASTA_FILTERED.md5_seq.sort | sort -k1,1V | awk '{printf( "%s\t%s\t%s\n", $1, $3, $5);}' | tee -a ${LOG_FILE}

  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


#
# Get sequence/chromosome sizes.
#
function make_chromosome_sizes_file()
{
  echo "Make chromosome sizes file..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  samtools faidx $FASTA_FILTERED
  REXP=`echo "$SEQUENCES_TO_KEEP_ANALYSIS" | sed 's/[ ][ ]*/|/g' | sed 's/^/(/' | sed 's/$/)/'`
  cat $FASTA_FILTERED.fai | awk '{if($1~/^'$REXP'$/) { printf( "%s\t%s\n",$1,$2); } }' > $CHROMOSOME_SIZES_FILE

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $CHROMOSOME_SIZES_FILE | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function make_whitelist_regions_file()
{
  echo "Make whitelist regions bed file..."  | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  cat $CHROMOSOME_SIZES_FILE \
  | awk 'BEGIN{OFS="\t"}{ print $1,"1",$2;}' > $WHITELIST_REGIONS_BED

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  cat $WHITELIST_REGIONS_BED | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}

