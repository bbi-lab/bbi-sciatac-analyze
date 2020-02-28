#!/bin/bash

function make_tss_file()
{
  echo "Make TSS bed file ${TSS_BED}.gz..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}

  Rscript $R_GENERATE_TSS_FILE $GTF_GZ ${TSS_BED}.temp
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count gene_biotype entries in ${TSS_BED}.temp by type..." | tee -a ${LOG_FILE}
  awk '{print$9}' ${TSS_BED}.temp | sort | uniq -c | sort -k 2,2 | awk '{printf( "%s\t%s\n", $1, $2 );}' > $TSS_FOUND_BIOTYPES_FILE
  cat $TSS_FOUND_BIOTYPES_FILE | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Expected gene_biotype entries in ${TSS_BED}.temp..."  | tee -a ${LOG_FILE}
  echo $SELECT_GENE_BIOTYPES | tr '|' '\n' | sort -k 1,1 > $TSS_SELECT_BIOTYPES_FILE
  cat $TSS_SELECT_BIOTYPES_FILE  | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Expected gene_biotypes not found in ${TSS_BED}.temp (may be either partial matches or the file uses different types)..." | tee -a ${LOG_FILE}
  join -1 1 -2 2 -v 1 $TSS_SELECT_BIOTYPES_FILE $TSS_FOUND_BIOTYPES_FILE | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "Filter TSS bed file ${TSS_BED}.temp entries..." | tee -a ${LOG_FILE}
  cat ${TSS_BED}.temp \
  | grep -i -E "${SELECT_GENE_BIOTYPES}" \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
  | sort -k1,1V -k2,2n -k3,3n \
  | bedtools intersect -a stdin -b $WHITELIST_REGIONS_BED \
  | uniq \
  | gzip > ${TSS_BED}.gz

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count entries by sequence name in ${TSS_BED}.gz" | tee -a ${LOG_FILE}
  zcat ${TSS_BED}.gz | awk '{print$1}' | sort | uniq -c | sort -k 2,2V | awk '{printf( "%s\t%s\n", $1, $2 );}' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo | tee -a ${LOG_FILE}
  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}


function make_gene_bodies_file()
{
  echo "Make gene bodies bed file ${GENE_BODIES_PLUS_UPSTREAM_BED}.gz..." | tee -a ${LOG_FILE}
  date '+%Y.%m.%d:%H.%M.%S' | tee -a ${LOG_FILE}

  Rscript $R_GENERATE_GENE_BODY_FILE $GTF_GZ ${GENE_BODIES_BED}.temp
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count gene_biotype entries in ${GENE_BODIES_BED}.temp by type..." | tee -a ${LOG_FILE}
  awk '{print$8}' ${GENE_BODIES_BED}.temp | sort | uniq -c | sort -k2,2 | awk '{printf( "%s\t%s\n", $1, $2 );}' > $GENE_BODIES_FOUND_BIOTYPES_FILE
  cat $GENE_BODIES_FOUND_BIOTYPES_FILE | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Expected gene_biotype entries in ${GENE_BODIES_BED}.temp..."  | tee -a ${LOG_FILE}
  echo $SELECT_GENE_BIOTYPES | tr '|' '\n' | sort -k 1,1 > $GENE_BODIES_SELECT_BIOTYPES_FILE
  cat $GENE_BODIES_SELECT_BIOTYPES_FILE  | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Expected gene_biotypes not found in ${GENE_BODIES_BED}.temp (may be either partial matches or the file uses different types)..." | tee -a ${LOG_FILE}
  join -1 1 -2 2 -v 1 $GENE_BODIES_SELECT_BIOTYPES_FILE $GENE_BODIES_FOUND_BIOTYPES_FILE | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}

  echo "Filter gene bodies bed file ${GENE_BODIES_BED}.temp entries..." | tee -a ${LOG_FILE}
  cat ${GENE_BODIES_BED}.temp \
  | grep -i -E "${SELECT_GENE_BIOTYPES}" \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
  | sort -k1,1V -k2,2n -k3,3n \
  | bedtools intersect -a stdin -b $WHITELIST_REGIONS_BED \
  | uniq \
  | gzip > ${GENE_BODIES_BED}.gz
  
  bedtools slop -i ${GENE_BODIES_BED}.gz -s -l 2000 -r 0 -g $CHROMOSOME_SIZES_FILE \
  | sort -k1,1V -k2,2n -k3,3n \
  | gzip > ${GENE_BODIES_PLUS_UPSTREAM_BED}.gz

  echo "CHECKPOINT" | tee -a ${LOG_FILE}
  echo "Count entries by sequence name in ${GENE_BODIES_BED}.gz" | tee -a ${LOG_FILE}
  zcat ${GENE_BODIES_BED}.gz | awk '{print$1}' | sort | uniq -c | sort -k 2,2V | awk '{printf( "%s\t%s\n", $1, $2 );}' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE} | tee -a ${LOG_FILE}


  echo | tee -a ${LOG_FILE}
  echo 'Done.' | tee -a ${LOG_FILE}
  echo | tee -a ${LOG_FILE}
}

