#!/bin/bash

function get_fasta_file
{
  rm -f $FASTA_GZ
  rm -f $FASTA
  rm -f ${CHECKSUMS}.dna
  rm -f ${README}.dna

  echo "Download and uncompress fasta file" | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  echo "URL is ${ENSEMBL_DNA_URL}" | tee -a ${LOG_FILE}
  echo "Fasta filename is $FASTA_GZ" | tee -a ${LOG_FILE}
  echo "Download fasta file..." | tee -a ${LOG_FILE}
  wget --no-verbose ${ENSEMBL_DNA_URL}/$FASTA_GZ 2>&1 | tee -a ${LOG_FILE}

  echo "Download ${CHECKSUMS}..." | tee -a ${LOG_FILE}
  wget --no-verbose ${ENSEMBL_DNA_URL}/${CHECKSUMS} 2>&1 | tee -a ${LOG_FILE}
  mv ${CHECKSUMS} ${CHECKSUMS}.dna 2>&1 | tee -a ${LOG_FILE}

  echo "Download ${README}..." | tee -a ${LOG_FILE}
  wget --no-verbose ${ENSEMBL_DNA_URL}/${README} 2>&1 | tee -a ${LOG_FILE}
  mv ${README} ${README}.dna 2>&1 | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo -n "Calculated $FASTA_GZ checksum is " | tee -a ${LOG_FILE}
  sum $FASTA_GZ | tee ${FASTA_GZ}.checksum | tee -a ${LOG_FILE}
  echo -n "Expected $FASTA_GZ checksum is " | tee -a ${LOG_FILE}
  grep $FASTA_GZ ${CHECKSUMS}.dna | awk '{print$1,$2}' 2>&1 | tee -a ${LOG_FILE}

  echo "Uncompress $FASTA file..." | tee -a ${LOG_FILE}
  zcat $FASTA_GZ > $FASTA

  echo | tee -a ${LOG_FILE}
  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


#
# Count the Ensembl sequence types in the fasta file.
#
function get_fasta_info()
{
  echo "Calculate $FASTA file sequence md5 checksums..." | tee -a ${LOG_FILE}
  $MD5_SEQ $FASTA > $FASTA.md5_seq

  echo "Get fasta sequence headers in $FASTA file..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}
  grep '^>' $FASTA > $FASTA.headers
 
  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Fasta sequence header in $FASTA file..." | tee -a ${LOG_FILE}
  cat $FASTA.headers | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Fasta sequence type counts in $FASTA file..." | tee -a ${LOG_FILE}
  awk '{print$4}' $FASTA.headers  | sort | uniq -c | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}

