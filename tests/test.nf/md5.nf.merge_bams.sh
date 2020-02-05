#!/bin/bash

lfil=`ls *.bam`

outf='md5.nf.merge_bams.out'
rm -f $outf

echo 'BEGIN_TABLE' | tee -a $outf
for fil in $lfil
do
  echo $fil | tee -a $outf
  samtools sort -n $fil | samtools view | md5sum | tee -a $outf
  echo | tee -a $outf
done
echo 'END_TABLE' | tee -a $outf
