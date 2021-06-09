#!/usr/bin/env python3

import sys
import re
import argparse


# timing information
# time with divert_mito_alignments.py
# s007> time filter.sh
#
# real    2m25.449s
# user    3m26.404s
# sys     0m18.139s
 
 
# time without divert_mito_alignments.py
# s007> time filter.sh
#
# real    1m57.977s
# user    2m28.572s
# sys     0m13.167s


# filter.sh script
# 
# in_file="bowtie2.sam"
# white_list="/net/bbi/vol1/data/genomes_stage/human/human_atac/whitelist_regions.bed"
# white_list_with_mt="/net/bbi/vol1/data/bge/genomes/genomes_stage.test.20210528/human/human_atac/whitelist_with_mt_regions.bed"
#  
# # cat $in_file | samtools view -h -L ${white_list_with_mt} -f3 -F12 -q10 | /net/gs/vol1/home/bge/git/bbi-sciatac-analyze/src/divert_mito_alignments.py -o bowtie_mito.sam | samtools view -L ${white_list} -f3 -F12 -q10 -bS - > bowtie2.whitelist.bam
# cat $in_file | samtools view -h -L ${white_list_with_mt} -f3 -F12 -q10 |  samtools view -L ${white_list} -f3 -F12 -q10 -bS - > bowtie2.whitelist.bam
# 

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to divert mitochondrial alignments to a SAM file.')
    parser.add_argument('-o', '--output_file', required=True, help='Name of SAM file in which to write mitochondrial alignments.')
    args = parser.parse_args()

    re1 = re.compile(r'@')
    re2 = re.compile(r'[!-?A-~]{1,254}\t[0-9]+\t(MT|Mt|HUMAN_MT|MOUSE_MT|MtDNA|mitochondrion_genome)')

    ofp = open(args.output_file, 'wt')

    for line in sys.stdin:
      if(re2.match(line)):
        ofp.write(line)
      elif(re1.match(line)):
        ofp.write(line)
        sys.stdout.write(line)
      else:
        sys.stdout.write(line)


