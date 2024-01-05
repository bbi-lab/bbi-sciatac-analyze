#!/bin/bash

infile="timeseriesPilot-aggregate.hash_reads.tsv"
sample_name="MLR_drugs"

# cat $infile | datamash -s -g 2,3,4 count 3 | datamash -g 1 sum 4 | awk -v S=$sample_name '{OFS="\t";i} {print S, $0}' > hashReads.per.cell.datamash.out

# cat $infile | sort | uniq | datamash -s -g 2,3,4 count 3 | datamash -g 1 sum 4 | awk -v S=$sample_name '{OFS="\t";i} {print S, $0}' > hashUMIs.per.cell.datamash.out

cat $infile | sort | uniq | datamash -s -g 1,2,4,5 count 3 > hashTable.datamash.out

