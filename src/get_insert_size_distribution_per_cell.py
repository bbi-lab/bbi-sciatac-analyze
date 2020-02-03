import pysam
import argparse
import collections
import io_functions as io

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Script to get insert size and total count per cell distributions per cell in long format.')
	parser.add_argument('fragments_file', help='BAM file with reads for each cell. Read names must be in the format "cellname:other_text" in order to allow cell IDs to be extracted.')
	parser.add_argument('output_file', help='Output file with insert size histograms for each cell.')
	parser.add_argument('--barcodes', help='Input file containing one cell ID per line. Only cell IDs included in this barcode list will be included in output files.')
	args = parser.parse_args()

	insert_sizes = {}

	# Read in the whitelist if provided
	barcodes = None
	if args.barcodes:
		barcodes = set([line.strip() for line in open(args.barcodes)])

	for line in io.open_file(args.fragments_file, 'rt'):
		chrom, start, end, cell, _ = line.strip().split('\t')

		# Get the insert size
		insert_size = int(end) - int(start)

		# Discount inferred insert sizes over 1000
		if insert_size >= 1000:
			continue

		# Track stats
		if barcodes and cell not in barcodes:
			continue

		if cell not in insert_sizes:
			insert_sizes[cell] = collections.Counter()

		insert_sizes[cell][insert_size] += 1

	# Output file with stats
	with open(args.output_file, 'w') as output_file:
		output_file.write('\t'.join(["cell", "insert_size", "fragment_count"]) + '\n')

		for cell in insert_sizes:
			counts = insert_sizes[cell]
			for insert_size in counts:
				output_file.write('\t'.join([cell, str(insert_size), str(counts[insert_size])]) + '\n')

