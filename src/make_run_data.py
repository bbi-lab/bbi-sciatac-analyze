#!/usr/bin/env python3

#
# Write run_data.js file for sci-ATAC experiment dashboard.
#

import sys
import argparse
import json
import csv
import re

import pprint


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


def get_sample_stats( args_json ):
    samples = get_samples( args_json )
    analyze_dir = args_json['analyze_dir']
    sample_stats = {}
    for sample in samples:
        sample_summary_file = "%s-called_cells_summary.stats.txt" % ( sample )
        read_tsv = csv.DictReader( open( sample_summary_file ), delimiter='\t' )
        row = next(read_tsv)
        stats = { 'Sample': row['sample'],
                  'Cell_threshold': row['cell_threshold'],
                  'Fraction_hs': '%.3f' % ( float( row['fraction_hs'] ) ),
                  'Fraction_tss': '%.3f' % ( float( row['fraction_tss'] ) ),
                  'Median_per_cell_frip': '%.3f' % ( float( row['median_per_cell_frip'] ) ),
                  'Median_per_cell_frit': '%.3f' % ( float( row['median_per_cell_frit'] ) ),
                  'Tss_enrichment': '%.3f' % ( float( row['tss_enrichment']) ),
                  'Sample_peaks_called': row['sample_peaks_called'],
                  'Total_merged_peaks': row['total_merged_peaks'],
                  'Total_reads': row['total_reads'],
                  'Fraction_reads_in_cells': '%.3f' % ( float( row['fraction_reads_in_cells'] ) ),
                  'Total_barcodes': row['total_barcodes'],
                  'Number_of_cells': row['number_of_cells'],
                  'Median_reads_per_cell': row['median_reads_per_cell'],
                  'Min_reads_per_cell': row['min_reads_per_cell'],
                  'Max_reads_per_cell': row['max_reads_per_cell'],
                  'Median_duplication_rate': '%.3f' % ( float( row['median_duplication_rate'] ) ),
                  'Median_fraction_molecules_observed': '%.3f' % ( float( row['median_fraction_molecules_observed'] ) ),
                  'Median_total_fragments': '%.3f' % ( float( row['median_total_fragments'] ) ),
                  'Total_deduplicated_reads': row['total_deduplicated_reads']
                }
        if( row.get( 'bloom_collision_rate' ) ):
            stats['Bloom_collision_rate'] = row['bloom_collision_rate']
        sample_stats.setdefault( row['sample'], stats )
    return( sample_stats )


def get_collision_rate( args_json ):
    samples = get_samples( args_json )
    analyze_dir = args_json['analyze_dir']
    sample_stats = {}
    for sample in samples:
        sample_summary_file = "%s-called_cells_summary.stats.txt" % ( sample )
        read_tsv = csv.DictReader( open( sample_summary_file ), delimiter='\t' )
        row = next(read_tsv)
        if( row.get( 'bloom_collision_rate' ) ):
            return( row['bloom_collision_rate'] )
    return( None )


def make_run_data_dict( args_json ):
    summarize_samples = get_sample_summarize( args_json )
    return( samples )


def write_run_data( out_file, run_name, samples, barnyard_collision_rate, sample_stats ):
    if( barnyard_collision_rate != None ):
      barn_collision = '%.1f%%' % ( float( barnyard_collision_rate ) * 100 )
    else:
      barn_collision = 'NA'
    run_data_dict = { 'run_name'    : run_name,
                      'cell_counts' : [ 0, 0 ],
                      'sample_list' : samples,
                      'barn_collision' : barn_collision,
                      'sample_stats'   : sample_stats }
    json_object = json.dumps( run_data_dict, indent = 2 )
    out_file.write( 'const run_data =\n' )
    out_file.write( json_object )
    return( 0 )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to write run_data.js file for sci-ATAC experiment dashboard.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to args.json input file in sci-ATAC analyze output directory.')
    parser.add_argument('-o', '--output_file', required=True, help='Name of output file.')
    args = parser.parse_args()
  
    args_json = read_args_json( open( args.input_file ) )
    run_name = get_run_name( args_json )
    samples = get_samples( args_json )
    sample_stats = get_sample_stats( args_json )
    collision_rate = get_collision_rate( args_json )
    write_run_data( open( args.output_file, 'w' ), run_name, samples, collision_rate, sample_stats )

