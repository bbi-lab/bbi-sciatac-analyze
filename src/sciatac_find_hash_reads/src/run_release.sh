#!/bin/bash

PGRM="/net/gs/vol1/home/bge/git/bbi-sciatac-analyze/src/sciatac_find_hash_reads/target/release/sciatac_find_hash_reads"
DATA_DIR="/net/bbi/vol1/data/regression_tests/sciATACseq/data_sources/fastqs/demux_fastqs/nextseq/sciplex-timeseriesPilot-AndrewMullen"


# echo "Run on uncompressed 10M read fastq files."
# time $PGRM -1 $DATA_DIR/timeseriesPilot-RUN001_L001_R1.10M.fastq -2 $DATA_DIR/timeseriesPilot-RUN001_L001_R2.10M.fastq -h $DATA_DIR/ATAC_HashBCList_Plate15.txt -o hits.out > /dev/null

echo "Run on compressed 10M read fastq files."
time $PGRM -1 $DATA_DIR/timeseriesPilot-RUN001_L001_R1.10M.fastq.gz -2 $DATA_DIR/timeseriesPilot-RUN001_L001_R2.10M.fastq.gz -h $DATA_DIR/ATAC_HashBCList_Plate15.txt -o hits.out > /dev/null

