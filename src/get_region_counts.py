from __future__ import print_function
import argparse
import sys
from io_functions import open_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to get count of reads within a given set of regions in the genome.')
    parser.add_argument('--transposition_sites_intersect', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for intersect.')
    parser.add_argument('--output_file', required=True, help='Output file of counts within region per cell.')
    args = parser.parse_args()

    # Read in an tabulate counts for each cell from piped in bed file
    cell_counts = {}
    for line in args.transposition_sites_intersect:
        cell = line.strip().split('\t')[3].split(':')[0]
            
        if cell not in cell_counts:
            cell_counts[cell] = 0
        cell_counts[cell] += 1

    # Make sure didn't get something empty
    if not cell_counts:
        raise ValueError('No counts found in specified regions for this sample: %s' % args.output_file)

    # Write out results
    with open(args.output_file, 'w') as output_file:
        output_file.write('\t'.join(['cell', 'count']) + '\n')

        for cell in sorted(cell_counts.keys()):
            output_file.write('\t'.join([cell, str(cell_counts[cell])]) + '\n')
