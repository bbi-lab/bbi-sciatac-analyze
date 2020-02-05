#!/bin/bash

lfil=`ls *-called_cells_stats.json *-called_cells.txt *-called_cells_whitelist.txt`

outf='md5.nf.call_cells.out'
rm -f $outf

echo '** unsorted' | tee -a $outf
for fil in $lfil
do
  echo $fil | tee -a $outf
  sort $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done

echo | tee -a $outf
echo | tee -a $outf

echo '** sorted' | tee -a $outf

echo 'BEGIN_TABLE' | tee -a $outf
for fil in $lfil
do
  echo $fil | tee -a $outf
  sort $fil | md5sum | tee -a $outf
  echo | tee -a $outf
done
echo 'END_TABLE' | tee -a $outf
