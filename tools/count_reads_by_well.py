#!/usr/bin/env python3

import argparse
import sys
import re


# TSV file with the contents
#
# cell    total   total_deduplicated  total_deduplicated_peaks    total_deduplicated_tss
# P2-B10_P3-A10_F08-rowD08-colF11 498 486 0   42
# P2-B06_P1-B10_G01-rowA01-colG08 1432    1382    4   77
# P4-G12_P2-D10_C08-rowB08-colC09 1322    1266    3   64

def read_count_report_file(filename):
  lig7 = {}
  lig5 = {}
  pcr7 = {}
  pcr5 = {}
  with open(filename, 'r') as fp:
    # Read header.
    nline = 0;
    line = fp.readline().rstrip();
    toks = line.split('\t')
    nline += 1
    # Check for header.
    if(toks[0] != 'cell' or toks[1] != 'total' or toks[2] != 'total_deduplicated'):
      print('Error: count_report file is missing header line.', file=sys.stderr)
      sys.exit(-1)
    for line in fp:
      nline += 1
      line = line.rstrip()
#      if(nline > 10):
#        break
      toks = line.split('\t')
      cell_name = toks[0]
      toks2 = cell_name.split('_')
      # Ligation 1 and 2.
      slig7 = toks2[0]
      slig5 = toks2[1]
      # PCR row and column
      toks3 = toks2[2].split('-')
      spcr7 = toks3[1]
      spcr5 = toks3[2]

      total_reads = int(toks[1])
      lig7[slig7] = lig7.setdefault(slig7, 0) + total_reads
      lig5[slig5] = lig5.setdefault(slig5, 0) + total_reads
      pcr7[spcr7] = pcr7.setdefault(spcr7, 0) + total_reads
      pcr5[spcr5] = pcr5.setdefault(spcr5, 0) + total_reads

  return( (lig7, lig5, pcr7, pcr5) )


if __name__ == '__main__':
  parser = argparse.ArgumentParser('Count reads by well from bbi-sciatac-analyze <sample_name>.count_report.txt file.')

  parser.add_argument('-i', '--input_file', required=True, help='Input filename.')
  args = parser.parse_args()

  (lig7, lig5, pcr7, pcr5) = read_count_report_file(args.input_file)

  print('barcode_type well read_count')
  for key in lig7.keys():
    print('lig7 %s %d' % (key, lig7[key]) )

  for key in lig5.keys():
    print('lig5 %s  %d' % (key, lig5[key]) )

  for key in pcr7.keys():
    print('pcr7 %s  %d' % (key, pcr7[key]) )

  for key in pcr5.keys():
    print('pcr5 %s  %d' % (key, pcr5[key]) )
  print()

