#!/bin/bash

DATA_DIR="/net/bbi/vol1/data/regression_tests/sciATACseq/data_sources/fastqs/demux_fastqs/nextseq/sciplex-timeseriesPilot-AndrewMullen"


echo "Run on uncompressed 1M read fastq files."
cargo run -- -1 $DATA_DIR/r1.fq.gz -2 $DATA_DIR/r2.fq.gz -h $DATA_DIR/ATAC_HashBCList_Plate15.txt -o hits.out

