#!/bin/bash

module purge
source /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load samtools/1.9
module load tbb/2019_U5 bowtie2/2.3.5.1
module load bedtools/2.28.0
module load R/3.6.1

#
# Executable paths.
#
MD5_SEQ="/net/bbi/vol1/data/src/sequtil/md5_seq"
FASTA_GETSEQS="/net/bbi/vol1/data/src/sequtil/fasta_getseqs"

#
# R script paths.
#
R_GENERATE_TSS_FILE="/net/bbi/vol1/data/bge/genomes/genomes_atac/ajh.scripts/sciatac_pipeline_data/generate_tss_file.R"
R_GENERATE_GENE_BODY_FILE="/net/bbi/vol1/data/bge/genomes/genomes_atac/ajh.scripts/sciatac_pipeline_data/generate_gene_body_file.R"

#
# Common file names.
#
FASTA=`echo $FASTA_GZ | sed 's/\.gz$//'`
FASTA_FILTERED=${FASTA}.filtered

CHECKSUMS="CHECKSUMS"
README="README"

FINAL_IDS_FILE="sequences_to_keep.txt"
CHROMOSOME_SIZES_FILE="chromosome_sizes.txt"

INDEX_PREFIX="$ORGANISM"

GTF=`echo $GTF_GZ | sed 's/\.gz$//'`

WHITELIST_REGIONS_BED="whitelist_regions.bed"
TSS_BED="tss.bed"
TSS_FOUND_BIOTYPES_FILE="tss.found_biotypes.lst"
TSS_SELECT_BIOTYPES_FILE="tss.select_biotypes.lst"

GENE_BODIES_BED="gene_bodies.bed"
GENE_BODIES_PLUS_UPSTREAM_BED="gene_bodies.plus_2kb_upstream.bed"
GENE_BODIES_FOUND_BIOTYPES_FILE="gene_bodies.found_biotypes.lst"
GENE_BODIES_SELECT_BIOTYPES_FILE="gene_bodies.select_biotypes.lst"
