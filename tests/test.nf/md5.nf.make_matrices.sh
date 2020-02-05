#!/bin/bash

outf='md5.nf.make_matrices.out'
rm -f $outf

lfil=`ls *gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *txt`
for fil in $lfil
do
  echo $fil | tee -a $outf
  md5sum $fil | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *genomic_windows.bed`
for fil in $lfil
do
  echo $fil | tee -a $outf
  md5sum $fil | tee -a $outf
  echo | tee -a $outf
done


echo | tee -a $outf
echo | tee -a $outf

echo '** sorted' | tee -a $outf

echo 'BEGIN_TABLE' | tee -a $outf
lfil=`ls *gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | sort | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *txt`
for fil in $lfil
do
  echo $fil | tee -a $outf
  sort $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *genomic_windows.bed`
for fil in $lfil
do
  echo $fil | tee -a $outf
  sort $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done
echo 'END_TABLE' | tee -a $outf
