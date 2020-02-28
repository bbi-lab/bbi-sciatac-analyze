#!/bin/bash

ORGANISM="macaque"

ENSEMBL_DNA_URL="ftp.ensembl.org:/pub/release-99/fasta/macaca_mulatta/dna"
FASTA_GZ="Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"

ENSEMBL_GTF_URL="ftp.ensembl.org:/pub/release-99/gtf/macaca_mulatta"
GTF_GZ="Macaca_mulatta.Mmul_10.99.gtf.gz"

#
# I see no sequence named MT in the macaque fasta file.
#
SEQUENCES_TO_KEEP_ALIGNER="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X Y MT"
SEQUENCES_TO_KEEP_ANALYSIS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X Y"
SELECT_GENE_BIOTYPES="protein|lncRNA|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV"

INDEX_PREFIX="$ORGANISM"
OUT_DIR="/net/bbi/vol1/data/bge/genomes/${ORGANISM}_atac"

LOG_FILE=$OUT_DIR/log.out


SCRIPT_DIR="/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze/scripts/genome"

source $SCRIPT_DIR/01.sciatac_genome.common_definitions.sh
source $SCRIPT_DIR/02.sciatac_genome.get_fasta_file.sh
source $SCRIPT_DIR/03.sciatac_genome.make_genome_files.sh
source $SCRIPT_DIR/04.sciatac_genome.make_aligner_indices.sh
source $SCRIPT_DIR/05.sciatac_genome.get_gtf_file.sh
source $SCRIPT_DIR/06.sciatac_genome.make_gene_bed_files.sh
source $SCRIPT_DIR/08.sciatac_genome.remove_unnecessary_files.sh

mkdir -p $OUT_DIR

pushd $OUT_DIR

#get_fasta_file
get_fasta_info
sequences_to_keep_named
filter_fasta_file
make_chromosome_sizes_file
make_whitelist_regions_file
make_aligner_index
get_gtf_file
get_gtf_info
make_tss_file
make_gene_bodies_file
#remove_unnecessary_files

popd

