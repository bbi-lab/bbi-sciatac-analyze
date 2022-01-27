#!/usr/bin/env python3

#
# Note: there is a Python module called PyPDF2, which is
#       supposed to merge PDF file. I am using the program
#       pdfunite because it is installed on the GS cluster.
#


#
# Merge summary PDF plot files.
#

import sys
import argparse
import json
import csv
import re
import os.path
import os


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]


def read_args_json( args_file ):
    """
    Read args.json file created in bbi-sciatac-analyze pipeline run.
    """
    args_json = json.load( args_file )
    return( args_json )


def dump_args_json( args_json ):
    """
    Diagnostic dump of args_json dictionary.
    """
    print('dump args json dictionary keys level 1')
    for dict_key in args_json.keys():
      print( 'dict_key: %s' % (dict_key))
  
    print('dump args json dictionary keys level 2 (runs)')
    for sample in args_json['samples']:
      print('sample: %s' % (sample))
    return( 0 )


def get_run_name( args_json ):
    analyze_dir = args_json['analyze_dir']
    parts = analyze_dir.split( '/' )
    return( parts[len(parts)-2] )


def get_samples( args_json ):
    samples = []
    for sample in args_json['samples']:
        samples.append( sample )
    return( samples )


def merge_summary_pdfs( output_filename, analyze_out, samples ):
    input_filename_list = []
    samples.sort(key=natural_sort_key)
    for sample in samples:
        filepath = '%s/%s/summarize_cell_calls/%s-called_cells_summary.pdf' % ( analyze_out, sample, sample )
        if( os.path.exists( filepath ) and os.access(filepath, os.R_OK ) ):
            input_filename_list.append( '%s/%s/summarize_cell_calls/%s-called_cells_summary.pdf' % ( analyze_out, sample, sample ) )
    command = 'pdfunite ' + ' '.join( input_filename_list ) + ' ' + output_filename
    os.system( command )
    return( 0 )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to merge summary PDF files.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to args.json input file in sci-ATAC analyze output directory.')
    parser.add_argument('-o', '--output_file', required=True, help='Name of output file.')
    args = parser.parse_args()
  
    args_json = read_args_json( open( args.input_file ) )
    samples = get_samples( args_json )
    
    merge_summary_pdfs( args.output_file, os.path.dirname( args.input_file ), samples )

