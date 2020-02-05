#!/bin/bash

outf='md5.ajh.tss_region_coverage.out'
rm -f $outf

echo '** unsorted' | tee -a $outf

lfil=`ls *tss_region_coverage.txt.gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

echo | tee -a $outf
echo | tee -a $outf

echo '** sorted' | tee -a $outf

echo 'BEGIN_TABLE' | tee -a $outf
lfil=`ls *tss_region_coverage.txt.gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | sort | md5sum | tee -a $outf
  echo | tee -a $outf
done
echo 'END_TABLE' | tee -a $outf


