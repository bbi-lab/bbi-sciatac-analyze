#!/usr/bin/env python3

import argparse
import sys

# Operation
#   o  read *-(mito.)duplicate_report.txt file, which has format
#        cell    total   total_deduplicated
#        P2-E12_P2-E12_F09-rowF09-colF06 19840   12748
#        P2-E10_P2-E10_D02-rowE02-colD05 353034  229138
#      Notes:
#        o  this is a tab-separated value file.
#        o  assume that the files may some cells that are not common
#           to both files
#        o  use argparse to get the name of each file and the output
#           file name
#        o  use dictionary to combine counts
#        o  the dictionary value is a list of counts
#             o  non-mito total count
#             o  non-mito total deduplicated count
#             o  mito total count
#             o  mito total deduplicated count
#             o  total count (?)
#             o  total deduplicated count (?)
#        o  used in bbi-sciatac-analyze main.nf process combineReadCountsProcess

def read_duplicate_report(combined_duplicate_report, file_name, offset):
    with open(file_name, 'r') as fp:
        tokens = fp.readline().split('\t')
        if(not all(list(tokens[0] == 'call', tokens[1] == 'total', tokens[2] == 'total_deduplicated'))):
            print('Error: unexpected header in file \'%s\'' % (file_name), file=sys.stderr)
            sys.exit(-1)
        for line in fp:
            tokens = line.split('\t')
            key = tokens[0]
            combined_duplicate_report.setdefault(key, [0, 0, 0, 0])
            combined_duplicate_report[key][0+offset] += int(tokens[1]
            combined_duplicate_report[key][1+offset] += int(tokens[2]


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to deduplicate position sorted sciATAC BAM file.')
    parser.add_argument('--input_duplicate_report', required=True, help='Input non-mitochondrial duplicate report file.')
    parser.add_argument('--input_mito_duplicate_report', required=True, help='Input mitochondrial duplicate report file..')
    parser.add_argument('--output_combined_duplicate_report', required=True, help='Output combined duplicate report file.')
    args = parser.parse_args()

    combined_duplicate_report = {}
    combined_duplicate_report = read_duplicate_report(combined_duplicate_report, args.input_duplicate_report, 0)
    combined_duplicate_report = read_duplicate_report(combined_duplicate_report, args.input_mito_duplicate_report, 2)
    write_duplicate_report(args.output_combined_duplicate_report, combined_duplicate_report)

