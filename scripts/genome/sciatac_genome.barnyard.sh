#!/bin/bash

ORGANISM="barnyard"

#
# Human definitions.
#
HUMAN_ENSEMBL_DNA_URL="ftp.ensembl.org:/pub/release-99/fasta/homo_sapiens/dna"
HUMAN_FASTA_GZ="Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

HUMAN_ENSEMBL_GTF_URL="ftp.ensembl.org:/pub/release-99/gtf/homo_sapiens"
HUMAN_GTF_GZ="Homo_sapiens.GRCh38.99.gtf.gz"

HUMAN_SEQUENCES_TO_KEEP_ALIGNER="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
HUMAN_SEQUENCES_TO_KEEP_ANALYSIS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
HUMAN_SELECT_GENE_BIOTYPES="protein|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV"

HUMAN_FINAL_IDS_FILE="human.sequences_to_keep.txt"
HUMAN_CHROMOSOME_SIZES_FILE="human.chromosome_sizes.txt"
HUMAN_WHITELIST_REGIONS_BED="human.whitelist_regions.bed"
HUMAN_TSS_BED="human.tss.bed"
HUMAN_TSS_FOUND_BIOTYPES_FILE="human.tss.found_biotypes.lst"
HUMAN_TSS_SELECT_BIOTYPES_FILE="human.tss.select_biotypes.lst"

#
# Mouse definitions.
#
MOUSE_ENSEMBL_DNA_URL="ftp.ensembl.org:/pub/release-99/fasta/mus_musculus/dna"
MOUSE_FASTA_GZ="Mus_musculus.GRCm38.dna.toplevel.fa.gz"

MOUSE_ENSEMBL_GTF_URL="ftp.ensembl.org:/pub/release-99/gtf/mus_musculus"
MOUSE_GTF_GZ="Mus_musculus.GRCm38.99.chr.gtf.gz"

MOUSE_SEQUENCES_TO_KEEP_ALIGNER="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT"
MOUSE_SEQUENCES_TO_KEEP_ANALYSIS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"
MOUSE_SELECT_GENE_BIOTYPES="protein|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV"

MOUSE_FINAL_IDS_FILE="mouse.sequences_to_keep.txt"
MOUSE_CHROMOSOME_SIZES_FILE="mouse.chromosome_sizes.txt"
MOUSE_WHITELIST_REGIONS_BED="mouse.whitelist_regions.bed"
MOUSE_TSS_BED="mouse.tss.bed"
MOUSE_TSS_FOUND_BIOTYPES_FILE="mouse.tss.found_biotypes.lst"
MOUSE_TSS_SELECT_BIOTYPES_FILE="mouse.tss.select_biotypes.lst"


#
# Common definitions.
#
OUT_DIR="/net/bbi/vol1/data/bge/genomes/${ORGANISM}_atac"

LOG_FILE=$OUT_DIR/log.out


SCRIPT_DIR="/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze/scripts/genome"

#
# Import functions.
#
source $SCRIPT_DIR/01.sciatac_genome.common_definitions.sh
source $SCRIPT_DIR/02.sciatac_genome.get_fasta_file.sh
source $SCRIPT_DIR/03.sciatac_genome.make_genome_files.sh
source $SCRIPT_DIR/04.sciatac_genome.make_aligner_indices.sh
source $SCRIPT_DIR/05.sciatac_genome.get_gtf_file.sh
source $SCRIPT_DIR/06.sciatac_genome.make_gene_bed_files.sh
source $SCRIPT_DIR/07.sciatac_genome.make_barnyard_files.sh
source $SCRIPT_DIR/08.sciatac_genome.remove_unnecessary_files.sh

#
# Begin processing.
#
mkdir -p $OUT_DIR

pushd $OUT_DIR

#
# Prepare human files.
#
ENSEMBL_DNA_URL=$HUMAN_ENSEMBL_DNA_URL
FASTA_GZ=$HUMAN_FASTA_GZ
FASTA=`echo $FASTA_GZ | sed 's/\.gz$//'`
FASTA_FILTERED=${FASTA}.filtered

ENSEMBL_GTF_URL=$HUMAN_ENSEMBL_GTF_URL
GTF_GZ=$HUMAN_GTF_GZ
GTF=`echo $GTF_GZ | sed 's/\.gz$//'`

SEQUENCES_TO_KEEP_ALIGNER=$HUMAN_SEQUENCES_TO_KEEP_ALIGNER
SEQUENCES_TO_KEEP_ANALYSIS=$HUMAN_SEQUENCES_TO_KEEP_ANALYSIS
SELECT_GENE_BIOTYPES=$HUMAN_SELECT_GENE_BIOTYPES

FINAL_IDS_FILE=$HUMAN_FINAL_IDS_FILE
CHROMOSOME_SIZES_FILE=$HUMAN_CHROMOSOME_SIZES_FILE
WHITELIST_REGIONS_BED=$HUMAN_WHITELIST_REGIONS_BED
TSS_BED=$HUMAN_TSS_BED
TSS_FOUND_BIOTYPES_FILE=$HUMAN_TSS_FOUND_BIOTYPES_FILE
TSS_SELECT_BIOTYPES_FILE=$HUMAN_TSS_SELECT_BIOTYPES_FILE

#
# Function calls for human.
#
get_fasta_file
get_fasta_info
get_gtf_file
get_gtf_info
sequences_to_keep_named
filter_fasta_file
make_chromosome_sizes_file
make_whitelist_regions_file
make_tss_file

#
# Rename files.
#
mv CHECKSUMS.dna human.CHECKSUMS.dna
mv README.dna human.README.dna
mv CHECKSUMS.gtf human.CHECKSUMS.gtf
mv README.dna human.README.gtf


#
# Prepare mouse files.
#
ENSEMBL_DNA_URL=$MOUSE_ENSEMBL_DNA_URL
FASTA_GZ=$MOUSE_FASTA_GZ
FASTA=`echo $FASTA_GZ | sed 's/\.gz$//'`
FASTA_FILTERED=${FASTA}.filtered

ENSEMBL_GTF_URL=$MOUSE_ENSEMBL_GTF_URL
GTF_GZ=$MOUSE_GTF_GZ
GTF=`echo $GTF_GZ | sed 's/\.gz$//'`

CHECKSUMS="mouse.CHECKSUMS"
README="mouse.README"

SEQUENCES_TO_KEEP_ALIGNER=$MOUSE_SEQUENCES_TO_KEEP_ALIGNER
SEQUENCES_TO_KEEP_ANALYSIS=$MOUSE_SEQUENCES_TO_KEEP_ANALYSIS
SELECT_GENE_BIOTYPES=$MOUSE_SELECT_GENE_BIOTYPES

FINAL_IDS_FILE=$MOUSE_FINAL_IDS_FILE
CHROMOSOME_SIZES_FILE=$MOUSE_CHROMOSOME_SIZES_FILE
WHITELIST_REGIONS_BED=$MOUSE_WHITELIST_REGIONS_BED
TSS_BED=$MOUSE_TSS_BED
TSS_FOUND_BIOTYPES_FILE=$MOUSE_TSS_FOUND_BIOTYPES_FILE
TSS_SELECT_BIOTYPES_FILE=$MOUSE_TSS_SELECT_BIOTYPES_FILE

#
# Function calls for mouse.
#

get_fasta_file
get_fasta_info
get_gtf_file
get_gtf_info
sequences_to_keep_named
filter_fasta_file
make_chromosome_sizes_file
make_whitelist_regions_file
make_tss_file

#
# Rename files.
#
mv CHECKSUMS.dna mouse.CHECKSUMS.dna
mv README.dna mouse.README.dna
mv CHECKSUMS.gtf mouse.CHECKSUMS.gtf
mv README.dna mouse.README.gtf


#
# Edit human and mouse files to change chromosome names.
#
barnyard_edit_fastas
barnyard_edit_chromosome_sizes
barnyard_edit_whitelist_regions
barnyard_edit_tss_bed

#
# Combine files and make aligner index.
#
barnyard_combine_files
#barnyard_make_aligner_index

#
# Clean up a bit.
#
#barnyard_remove_unnecessary_files

popd

