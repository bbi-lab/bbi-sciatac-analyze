#!/bin/bash

# cat $WORKING_DIR/parse_hash_out/Parsed_Hash_Table \
#     | $DATAMASH_PATH -s -g 2,3,4 count 3 \
#     | $DATAMASH_PATH -g 1 sum 4 \
#     | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
#     > $WORKING_DIR/hashRDS/hashReads.per.cell
# 
# cat $WORKING_DIR/parse_hash_out/Parsed_Hash_Table \
#     | uniq \
#     | $DATAMASH_PATH -s -g 2,3,4 count 3 \
#     | $DATAMASH_PATH -g 1 sum 4 \
#     | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
#     > $WORKING_DIR/hashRDS/hashUMIs.per.cell
# 
# Rscript $SCRIPTS_DIR/knee-plot.R            \
#     $WORKING_DIR/hashRDS/hashUMIs.per.cell          \
#     $WORKING_DIR/hashRDS
# 
# cat $WORKING_DIR/parse_hash_out/Parsed_Hash_Table \
#     | uniq \
#     | $DATAMASH_PATH -s -g 1,2,4,5 count 3  \
#     > $WORKING_DIR/hashRDS/hashTable.out
# 
# paste $WORKING_DIR/hashRDS/hashUMIs.per.cell $WORKING_DIR/hashRDS/hashReads.per.cell\
#     | cut -f 1,2,6,3 \
#     | awk 'BEGIN {OFS="\t";} {dup = 1-($3/$4); print $1,$2,$3,$4,dup;}' \
#     > $WORKING_DIR/hashRDS/hashDupRate.txt


infile="timeseriesPilot-aggregate.hash_reads.tsv"
sample_name="MLR_drugs"

cat $infile | datamash -s -g 2,3,4 count 3 | datamash -g 1 sum 4 | awk -v S=$sample_name '{OFS="\t";} {print S, $0}' > hashReads.per.cell.datamash.out

cat $infile | sort | uniq | datamash -s -g 2,3,4 count 3 | datamash -g 1 sum 4 | awk -v S=$sample_name '{OFS="\t";} {print S, $0}' > hashUMIs.per.cell.datamash.out

cat $infile | sort | uniq | datamash -s -g 1,2,4,5 count 3 > hashTable.datamash.out

