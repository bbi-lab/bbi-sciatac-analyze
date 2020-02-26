#!/bin/bash

function remove_unnecessary_files()
{
  rm -f $FASTA_GZ $FASTA ${FASTA}.filtered
  rm -f $GTF_GZ
  rm -f ${TSS_BED}.temp
  rm -r ${GENE_BODIES_BED}.temp
}
