#!/usr/bin/env python3

import matplotlib
matplotlib.use('pdf')
import scrublet as scr
import scipy.io
import scipy.sparse
import numpy
import numpy.ma
from PIL import Image, ImageDraw, ImageFont
import os
import sys
import re
import warnings
import traceback
import argparse

#
# Notes:
#   o  apply umi_cutoff in filter_counts_matrix() that is consistent with the
#      umi_cutoff value used in reduce_dimensions.R in order to produce consistent
#      numbers of cells (columns).
#   o  I have the impression that scrublet is vulnerable to internal
#      resulting from things like division by zero and sqrt of x < 0.
#   o  it appears that python issues a runtime warning for some of
#      these errors rather than stopping so I am raising an exceptioni
#      in these cases
#

def handle_warning(message, category, filename, lineno, file=None, line=None):
    print( 'Scrublet: stop on warning \'%s\' in %s at line %s' % ( message, filename, lineno ), file=sys.stderr )
    raise ValueError(message)
    return( 0 )


def read_col_file(col_file):
  col = []
  with open(col_file, 'r') as fp:
    for line in fp:
      col.append(line.rstrip())
  return(col)


def filter_counts_matrix(mat_in, outlier_filter, umi_cutoff, col_names_file):
    print('run_scrublet.py: filter_counts_matrix: begin')
    # start with COO format matrix
    if(not scipy.sparse.isspmatrix_coo(mat_in)):
        mat_in = mat_in.tocoo()
    # read column names
    col_in = read_col_file(col_names_file)
    # binarize matrix using directly the (non-zero) m.data attribute
    # (from snapATAC)
    cutoff = numpy.percentile(a=mat_in.data, q=100 - outlier_filter, axis=None)
    mat_in.data[mat_in.data > cutoff] = 0
    mat_in.data[mat_in.data > 1] = 1
    # find cells with no more than umi_cutoff counts
    colsum = mat_in.sum(0)[0, :]
    keep_col = colsum > umi_cutoff
    # subset matrix and column names 
    mat_out = mat_in.tocsc()[:, keep_col.tolist()[0]]
    col_out = [i for (i, v) in zip(col_in, keep_col.tolist()[0]) if v]
    print('run_scrublet.py: filter_counts_matrix: end')
    return(mat_out, col_out)


def run_scrublet(sample_name, counts_matrix):
    print('run_scrublet.py: run_scrublet: begin')
    warnings.showwarning = handle_warning

    if(numpy.size(counts_matrix, 0) == 0 or numpy.size(counts_matrix, 1) == 0):
        filename = args.sample_name + "-scrublet_hist.png"
        image = Image.new(mode = "RGB", size = (800,600), color = "white")
        draw = ImageDraw.Draw(image) 
        draw.text((50,50), "Scrublet failed. This is generally because there aren't enough cells with sufficient reads.\n", fill = "black")
        return(-1)

    if(not scipy.sparse.isspmatrix_csc(counts_matrix)):
        counts_matrix = counts_matrix.T.tocsc()
    else:
        counts_matrix = counts_matrix.T

    # count_matrix
    #   rows: cells
    #   cols: genes
    scrub = scr.Scrublet(counts_matrix)

    try:
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        scrub.plot_histogram()[0].savefig(args.sample_name + "-scrublet_hist.png")
        all_scores = numpy.vstack((doublet_scores, predicted_doublets))
        all_scores = numpy.transpose(all_scores)
        numpy.savetxt(args.sample_name + "-scrublet_table.csv", all_scores, delimiter=",", fmt='%.8e,%d')
    except (ZeroDivisionError, FloatingPointError, ValueError) as eobj:
        tb_str = traceback.format_exc()
        print('%s' % ( tb_str ), file=sys.stderr)
        temp = numpy.array(["NA"] * numpy.size(counts_matrix, 0))
        all_scores = numpy.vstack((temp, temp))
        all_scores = numpy.transpose(all_scores)
        filename = args.sample_name + "-scrublet_hist.png"
        image = Image.new(mode = "RGB", size = (800,600), color = "white")
        draw = ImageDraw.Draw(image)
        draw.text((50,50), "Scrublet failed. This is generally because there aren't enough cells with sufficient reads.\n\nFailure message:\n\n" + tb_str, fill = "black")
        image.save(filename)
        numpy.savetxt(args.sample_name + "-scrublet_table.csv", all_scores, fmt="%s", delimiter=",")
    except (AttributeError) as eobj:
        tb_str = traceback.format_exc()
        print('%s' % ( tb_str ), file=sys.stderr)
        predicted_doublets = scrub.call_doublets(threshold=0.15)
        scrub.plot_histogram()[0].savefig(args.sample_name + "-scrublet_hist.png")
        all_scores = numpy.vstack((doublet_scores, predicted_doublets))
        all_scores = numpy.transpose(all_scores)
        numpy.savetxt(args.sample_name + "-scrublet_table.csv", all_scores, delimiter=",", header='doublet_score,doublet')
    print('run_scrublet.py: run_scrublet: end')
    return( 0 )


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Run scrublet.')

    parser.add_argument('--sample_name', required=True, help='Sample name (for naming).')
    parser.add_argument('--mat_file', required=True, help='input matrix file name.')
    parser.add_argument('--umi_cutoff', type=int, required=True, help='umi filter cutoff.')
    args = parser.parse_args()

    col_names_file = re.sub('[.]gz$', '', args.mat_file)
    col_names_file = re.sub('[.]mtx$', '', col_names_file) + '.columns.txt'

    counts_matrix = scipy.io.mmread(args.mat_file)
    # note: numpy percentile takes percentile value whereas R quantile takes 'quantile' value: differ by factor of 100
    outlier_filter = 1.0e-1
    counts_matrix, col_names = filter_counts_matrix(counts_matrix, outlier_filter, args.umi_cutoff, col_names_file)
    run_scrublet(args.sample_name, counts_matrix)

    peak_matrix_cols_out_file = '%s-scrublet_columns.txt' % (args.sample_name)
    with open(peak_matrix_cols_out_file, 'w' ) as fp:
        fp.write('\n'.join(col_names) + '\n')

