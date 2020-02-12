#!/bin/bash

#
# What have we got here?
#
# hg19 json file entry
#        "name": "hg19",
#         "tss": "/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/hg19/tss.bed.gz",
#         "gene_score_bed": "/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/hg19/gene_bodies.plus_2kb_upstream.bed.gz",
#         "macs_genome": "hs",
#         "bowtie_index": "/net/shendure/vol10/projects/scATAC/nobackup/genomes/hg19/bowtie/hs37d5",
#         "whitelist_regions": "/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/hg19/whitelist_regions.bed",
#         "chromosome_sizes": "/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/hg19/chromosome_sizes.txt",
#         "motifs": "/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm"
# ,
#         "fasta": "/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/hg19/reference/hg19.fa"
#
# Notes:
#   o  macs_genome is the effective genome length: apparently some values are built-in
#

########################################
# Initialize variables
########################################

# Build-specific files
GENOME_VERSION="hg38"
GTF_FILE=""
SPECIES_NAME=""
SPECIES_COMMON_NAME=""

GENOME_DIRECTORY="/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze/genomes/scripts/genomes/hsapiens"
BOWTIE2_INDEX_DIRECTORY="$GENOME_DIRECTORY/genome_bowtie"


#
# Gencode full assembly URLs.
# REFERENCE_FASTA_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.p13.genome.fa.gz"
# REFERENCE_FASTA_GZ="GRCh38.p13.genome.fa.gz"
# REFERENCE_FASTA="GRCh38.p13.genome.fa"

# REFERENCE_GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz"
# REFERENCE_GTF_GZ="gencode.v33.annotation.gtf.gz"
# REFERENCE_GTF="gencode.v33.annotation.gtf"

# Gencode primary assembly URLs.
REFERENCE_FASTA_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz"
REFERENCE_FASTA_GZ="GRCh38.primary_assembly.genome.fa.gz"
REFERENCE_FASTA="GRCh38.primary_assembly.genome.fa"

REFERENCE_GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.primary_assembly.annotation.gtf.gz"
REFERENCE_GTF_GZ="gencode.v33.primary_assembly.annotation.gtf.gz"
REFERENCE_GTF="gencode.v33.primary_assembly.annotation.gtf"


########################################
# Define genome-specific variables.
########################################

# Script directory.
SCRIPT_DIR="/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze/genomes/scripts"


########################################
# Load required modules.
########################################
module purge
source /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load samtools/1.9
module load bowtie2/2.2.3


########################################
# Make required directories.
########################################
mkdir -p ${GENOME_DIRECTORY}/fasta ${GENOME_DIRECTORY}/gtf ${GENOME_DIRECTORY}/pfms


########################################
# Download required files.
########################################
echo 'Downloading required files...'

# reference fasta file.
if [ ! -f "$GENOME_DIRECTORY/fasta/REFERENCE_FASTA_GZ" ]
then
  pushd $GENOME_DIRECTORY/fasta
  wget $REFERENCE_FASTA_URL
  popd
fi

# reference annotation GTF file.
if [ ! -f "$GENOME_DIRECTORY/fasta/REFERENCE_FASTA_GZ" ]
then
  pushd $GENOME_DIRECTORY/gtf
  wget REFERENCE_GTF_URL
  popd
fi

# pfms
pushd $GENOME_DIRECTORY/pfms
wget http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt
mv JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm
popd

# record file sources to file for distribution with sample data
# xxx


########################################
# Genome files.
########################################
echo 'Prepping genome files for ${GENOME_VERSION}...'
pushd ${GENOME_DIRECTORY}/fasta

# REFERENCE (FASTA)
# sequence lengths 
# Notes:
#   o  filter out unwanted sequences
gunzip -c $REFERENCE_FASTA_GZ > $REFERENCE_FASTA
samtools faidx $REFERENCE_FASTA
awk 'BEGIN{FS="\t"}{printf( "%s\t%s\n",$1,$2 )}' $REFERENCE_FILE \
    | grep -v chrM \
    | grep -v random \
    > chromosome_sizes.txt

if [ ! -f reference/hg19.fa ]; then
    echo '   reference...'
    mkdir -p reference
    zcat ~ajh24/common_data/ucsc/goldenPath/hg19/bigZips/hg19.fa.gz > reference/hg19.fa
    # Also index the fasta
    python -c "import pyfasta; test = pyfasta.Fasta('reference/hg19.fa')"
fi


# BOWTIE2 REFERENCE
# make bowtie2 reference file
PARS=""
bowtie2-build $PARS $REFERENCE_LIST $BOWTIE2_INDEX_DIRECTORY

# WHITELIST REGIONS (All standard chromosomes from end to end)
if [ ! -f whitelist_regions.bed ]; then
    echo '   whitelist regions...'
    cat ~ajh24/common_data/ucsc/goldenPath/hg19/bigZips/hg19.chrom.sizes \
    | grep -v chrM \
    | grep -v random \
    > chromosome_sizes.txt

    cat chromosome_sizes.txt \
    | awk 'BEGIN{OFS="\t"}{ print $1,"1",$2;}' > whitelist_regions.bed


fi


# TSS AND GENE BODY + 2KB DEFINITIONS
if [ ! -f tss.bed.gz ]; then
    echo '   gtf...'
    mkdir -p gtf
    curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz > gtf/gencode.v29lift37.annotation.gtf.gz

    Rscript ~ajh24/common_data/sciatac_pipeline_data/generate_tss_file.R gtf/gencode.v29lift37.annotation.gtf.gz tss.temp.bed
    Rscript ~ajh24/common_data/sciatac_pipeline_data/generate_gene_body_file.R gtf/gencode.v29lift37.annotation.gtf.gz gene_bodies.temp.bed

    cat tss.temp.bed \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b whitelist_regions.bed \
    | uniq \
    | gzip > tss.bed.gz

    cat gene_bodies.temp.bed \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b whitelist_regions.bed \
    | uniq \
    | gzip > gene_bodies.bed.gz

    bedtools slop -i gene_bodies.bed.gz -s -l 2000 -r 0 -g chromosome_sizes.txt \
    | sort -k1,1V -k2,2n -k3,3n \
    | gzip > gene_bodies.plus_2kb_upstream.bed.gz

    rm tss.temp.bed
    rm gene_bodies.temp.bed
fi

popd

echo 'Done.'


