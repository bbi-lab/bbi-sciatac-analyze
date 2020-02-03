
from __future__ import print_function
import argparse
import os
import gzip
import subprocess
import itertools
import collections
import io_functions
import sys
import scipy.sparse as sp


def open_file(file_name, mode=None):
    if '.gz' == file_name[-3:]:
        if mode is None:
            mode = 'rt'
        return gzip.open(file_name, mode)
    else:
        if mode is None:
            mode = 'r'
        return open(file_name, mode)

def load_cell_whitelist(file_object):
    whitelist = dict()
    for i,line in enumerate(file_object):
         whitelist[line.strip()] = i

    return whitelist

def load_feature_whitelist(file_object):
    bed_entries = None
    feature_whitelist = dict()

    i = 0
    for line in file_object:
        entries = line.strip().split('\t')
        bed_entries = len(entries)

        if bed_entries == 3:
            if tuple(entries) not in feature_whitelist:
                feature_whitelist[tuple(entries)] = i
                i += 1
        elif bed_entries > 3:
            if entries[3] not in feature_whitelist:
                feature_whitelist[entries[3]] = i
                i += 1
        else:
            raise ValueError('Intervals must be BED file with at least 3 columns.')
    return feature_whitelist, bed_entries

def get_feature_cell_counts(file_object, feature_keys, cell_keys, interval_columns):
    peak_cell_counts = collections.Counter()

    # Group by each interval or each interval + group if there are four or more columns
    key_columns = min(4, interval_columns)

    for _, group in itertools.groupby(file_object, key=lambda line: tuple(line.strip().split('\t')[0:key_columns])):
        for line in list(group):
            entries = line.strip().split('\t')
            cell = entries[interval_columns + 3].split(":")[0]

            if cell in cell_keys:
                feature = None
                if interval_columns == 3:
                    feature = tuple(entries[0:3])
                else:
                    feature = entries[3]

                peak_cell_counts[(feature_keys[feature], cell_keys[cell])] += 1
    
    feature_indices = []
    cell_indices = []
    values = []
    for (feature, cell),v in peak_cell_counts.items():
        feature_indices.append(feature)
        cell_indices.append(cell)
        values.append(v)

    return feature_indices, cell_indices, values

def get_value_sorted_keys(d):
    return sorted(d.keys(), key=lambda x: d[x])


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to take BAM and BED intervals and generate a sparse peak by cell matrix. Note that if the BED file of intervals has 4 or more columns, the fourth column is used as the feature and the score will be the aggregate of all intervals belonging to that group.')
    parser.add_argument('--transposition_sites_intersect', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for intersect of transposition_sites_bed with features from which you are making matrix (e.g. bedtools intersect -sorted -a intervals.bed -b transposition_sites.bed -wa -wb).')
    parser.add_argument('--intervals', required=True, help='Sci-ATAC BAM file. Intervals will be features if three columns. If four columns, the fourth column will be used as the feature and score will be aggregate of all intervals belonging to that group.')
    parser.add_argument('--cell_whitelist', required=True, help='File with list of cell IDs allowed to output to final matrix.')
    parser.add_argument('--matrix_output', help='.mtx or .mtx.gz formatted output file.')
    args = parser.parse_args()

    # Load whitelists
    cell_whitelist = load_cell_whitelist(open_file(args.cell_whitelist))
    features, interval_columns = load_feature_whitelist(open_file(args.intervals))

    print('Building matrix...')
    feature_indices, cell_indices, values = get_feature_cell_counts(args.transposition_sites_intersect, features, cell_whitelist, interval_columns)
    feature_names = get_value_sorted_keys(features)
    if interval_columns == 3:
        feature_names = ['_'.join(feature) for feature in feature_names]
    else:
        feature_names = feature_names

    cell_names = get_value_sorted_keys(cell_whitelist)

    print('Writing output...')
    sparse_mat = sp.csc_matrix((values, (feature_indices, cell_indices)), shape=(len(feature_names), len(cell_names)), dtype='int32')
    io_functions.write_mtx_file(sparse_mat, feature_names, cell_names, args.matrix_output)
    print("Done.")