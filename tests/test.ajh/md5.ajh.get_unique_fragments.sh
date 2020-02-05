#!/bin/bash

outf='md5.ajh.get_unique_fragments.out'
rm -f $outf

echo '** unsorted' | tee -a $outf

lfil=`ls *bed.gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *duplicate_report.txt`
for fil in $lfil
do
  echo $fil | tee -a $outf
  md5sum $fil | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *fragments.txt.gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *insert_sizes.txt`
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
lfil=`ls *bed.gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | sort | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *duplicate_report.txt`
for fil in $lfil
do
  echo $fil | tee -a $outf
  sort $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *fragments.txt.gz`
for fil in $lfil
do
  echo $fil | tee -a $outf
  zcat $fil | sort | md5sum | tee -a $outf
  echo | tee -a $outf
done

lfil=`ls *insert_sizes.txt`
for fil in $lfil
do
  echo $fil | tee -a $outf
  sort $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done
echo 'END_TABLE' | tee -a $outf

