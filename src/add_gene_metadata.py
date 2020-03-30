#!/usr/bin/env python3

import argparse
import sys


def open_file( file_name, mode = None ):
    if '.gz' == file_name[-3:]:
        if mode is None:
            mode = 'rt'
        return gzip.open(file_name, mode)
    else:
        if mode is None:
            mode = 'r'
        return open(file_name, mode)


def load_row_names(file_object):
    row_names = []
    for line in file_object:
        fields = line.strip().split('\t')
        gene_id = fields[0]
        row_names.append(fields[0])
    return(row_names)


def load_gene_metadata(file_object):
    gene_metadata = {}
    for line in file_object:
        fields = line.strip().split('\t')
        gene_id = fields[0]
        gene_name = fields[1]
        gene_biotype = fields[2]
        gene_metadata.setdefault(gene_id, {})
        gene_metadata[gene_id]['gene_name'] = gene_name
        gene_metadata[gene_id]['gene_biotype'] = gene_biotype
    return(gene_metadata)


if __name__ == '__main__':
    parser = argparse.ArgumentParser( 'Script to add meta data columns to promoter matrix rows file.' )
    parser.add_argument('--in_promoter_matrix_row_name_file', required=True, help='Input file with promoter matrix official gene names.')
    parser.add_argument('--gene_metadata_file', required=True, help='Input file with gene metadata. The first column has the official gene names. Additional rows have metadata such as short gene names and gene biotypes.')
    parser.add_argument('--out_promoter_matrix_row_name_file', required=True, help='Output file with promoter matrix official gene names and metadata columns.')
    args = parser.parse_args()

    print('in_promoter_matrix_row_name_file: %s' % (args.in_promoter_matrix_row_name_file))
    print('gene_metadata_file: %s' % (args.gene_metadata_file))
    
    gene_metadata = load_gene_metadata(open_file(args.gene_metadata_file))
    row_names = load_row_names(open_file(args.in_promoter_matrix_row_name_file))
    out_file = open(args.out_promoter_matrix_row_name_file,'w')
    for row in row_names:
        if( row in gene_metadata ):
            print( '%s\t%s\t%s' % (row, gene_metadata[row]['gene_name'], gene_metadata[row]['gene_biotype']), file=out_file)
        else:
            print( '%s\t.\t.' % (row), file=out_file)
