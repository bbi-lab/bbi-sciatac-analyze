from __future__ import print_function
import pysam
import sys
import argparse
from collections import OrderedDict

MAX_INSERT_SIZE = 1000

def write_bam_bed(fragments_dict, file_handle, chromosome, flanking_distance=100):
    """
    Write out a set of intervals surrounding each transposition site accounted for by fragments in fragments dict to BED format.
    One line per transposition site.

    Args:
        fragments_dict (OrderedDict): dictionary of (cell, start, end) for fragments to the number of times that fragment was seen
        file_handle (file handle): open file handle to output file
        input_bam (Samfile): handle to BAM file from pysam
        reference_id (int): reference id as specified by pysam's tid property for AlignedSegments
        flanking_distance (int): output regions will be centered on transposition site and extended in either direction by this many bp
    """

    # Get all endpoints from all fragments
    endpoints = []
    for fragment in fragments_dict:
        cell, start, stop = fragment
        endpoints.append((start, cell))
        endpoints.append((stop, cell))

    # NOTE this is necessary because BAMs are not sorted by position of their transposition site,
    # transposition site != the position property of the read itself, which is what BAM is sorted by.
    # Could also just write out unsorted and sort later, instead just sorting here for now.
    endpoints.sort(key=lambda x: x[0])

    for coord, cell in endpoints:
        start = coord - flanking_distance
        end = coord + flanking_distance
        file_handle.write('\t'.join([chromosome, str(max(0, start)), str(end), cell]) + '\n')

def write_fragments(fragments_dict, file_handle, chromosome):
    """
    Write current set of fragments out to file. Fragments dict must be OrderedDict such that elements were inserted in position order.

    Args:
        fragments_dict (OrderedDict): dictionary of (cell, start, end) for fragments to the number of times that fragment was seen
        file_handle (file handle): open file handle to output file
        chromosome (string): chromosome for all entries (only one chromosome supported at a time)
    """
    for fragment in fragments_dict:
        cell, start, stop = fragment
        duplicate_count = observed_fragments[fragment]
        file_handle.write('\t'.join([chromosome, str(start), str(stop), cell, str(duplicate_count)]) + '\n')

def update_cell_duplicate_counts(cell_duplicate_counts, fragments_dict):
    """
    Updates running dict of cell duplicate counts using fragments dict.

    Args:
        fragments_dict (OrderedDict): dictionary of (cell, start, end) for fragments to the number of times that fragment was seen
        cell_duplicate_counts (dict): dict mapping cell ID to [deduplicated_counts, total_counts]
    
    Returns:
        dict: updated cell_duplicate_counts accounting for entries in fragments_dict
    """
    for fragment in fragments_dict:
        duplicate_counts = fragments_dict[fragment]
        cell = fragment[0]

        if cell not in cell_duplicate_counts:
            cell_duplicate_counts[cell] = [0, 0]
        
        cell_duplicate_counts[cell][0] += 1
        cell_duplicate_counts[cell][1] += duplicate_counts

    return cell_duplicate_counts

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to deduplicate position sorted sciATAC BAM file.')
    parser.add_argument('input_bam', help='Input BAM file to be deduplicated.')
    parser.add_argument('--output_bam', required=False, help='Output BAM file, deduplicated.')
    parser.add_argument('--fragments', required=False, help='Output file replicating the 10x fragments.tsv.gz format for visualization purposes, etc.')
    parser.add_argument('--transposition_sites_bed', required=False, help='Output file for peak calling and matrix generation. BED file with each region centered around a transposition site annotated with cell ID.')
    parser.add_argument('--duplicate_read_counts', required=False, help='Output file specifying the number of reads per cell (not fragments) in duplicated and deduplicated BAM files.')
    parser.add_argument('--insert_sizes', required=False, help='Output file with insert size distribution for the sample.')
    parser.add_argument('--flanking_distance', default=1, type=int, help='Optional and only used when specifying file for transposition_sites_out. Will make transposition site regions N bp wide centered on the site itself.')
    args = parser.parse_args()

    readsin = pysam.AlignmentFile(args.input_bam, "rb")

    if not (args.output_bam or args.fragments or args.transposition_sites_bed or args.duplicate_read_counts):
        raise ValueError("One or more of --output_bam, --fragments_out, --transposition_sites_bed, or --duplicate_read_counts should be specified.")

    if args.output_bam:
        readsout = pysam.AlignmentFile(args.output_bam, "wb", template=readsin)
    if args.fragments:
        fragments_out = open(args.fragments, 'w')
    if args.transposition_sites_bed:
        transposition_sites_out = open(args.transposition_sites_bed, 'w')

    last_reference = None
    observed_fragments = OrderedDict() # track duplicate counts
    read_names_written = set() # helps maintain writing matched read pairs when writing BAM
    cell_duplicate_counts = dict()
    insert_sizes = [0] * (MAX_INSERT_SIZE + 1)
    chromosome_name = None

    for read in readsin:
        reference_id = read.reference_id

        if chromosome_name is None:
            chromosome_name = chromosome = readsin.get_reference_name(reference_id)
    
        # If hit new reference ID, write per-chromosome fragment data out and reset
        if last_reference != reference_id and last_reference is not None:
            if args.fragments:
                write_fragments(observed_fragments, fragments_out, chromosome_name)
            if args.transposition_sites_bed:
                write_bam_bed(observed_fragments, transposition_sites_out, chromosome_name, flanking_distance=args.flanking_distance)

            chromosome_name = readsin.get_reference_name(reference_id)
            cell_duplicate_counts = update_cell_duplicate_counts(cell_duplicate_counts, observed_fragments)
            observed_fragments = OrderedDict()
            read_names_written = set()
        
        last_reference = reference_id

        # Now process read (note BowTie2 should not soft-clip, which allows direct use of pysam)
        # have compared this to compute_five_prime_coords in 10x pipeline just to be sure and always equal
        # for our BAM files
        cell = read.qname.split(':')[0]

        if read.tlen < 0:
            fragment_start = read.mpos
            fragment_end = read.mpos - read.tlen
        else:
            fragment_start = read.pos
            fragment_end = read.pos + read.tlen
    
        # Given extracted info, this ID should identify unique molecules with same endpoints per-chromosome
        fragment_id = (cell, fragment_start, fragment_end)

        if fragment_id not in observed_fragments:
            # first time fragment seen, track fragment/specific readname and write out
            observed_fragments[fragment_id] = 0
            
            # Track insert sizes if requested
            if args.insert_sizes:
                insert_size = fragment_end - fragment_start
                if insert_size <= MAX_INSERT_SIZE:
                    insert_sizes[insert_size] += 1

            if args.output_bam:
                read_names_written.add(read.qname)
                readsout.write(read)
        else:
            if args.output_bam and read.qname in read_names_written:
                # Mate of this read was written, so write it out specifically to maintain proper pairing
                readsout.write(read)

        # Track the number of duplicate fragments, but don't double count
        # since this tracks fragments, not reads
        if read.is_read1:
            observed_fragments[fragment_id] += 1
                    

    # Write out any more fragments as needed
    if observed_fragments:
        print('Finished %s and writing fragments to file...' % reference_id)
        if args.fragments:
            write_fragments(observed_fragments, fragments_out, chromosome_name)
        if args.transposition_sites_bed:
            write_bam_bed(observed_fragments, transposition_sites_out, chromosome_name, flanking_distance=args.flanking_distance)

        cell_duplicate_counts = update_cell_duplicate_counts(cell_duplicate_counts, observed_fragments)

    # Close out the output files
    if args.fragments:
        fragments_out.close()
    if args.transposition_sites_bed:
        transposition_sites_out.close()

    readsin.close()

    # Check to make sure something was read
    if not cell_duplicate_counts:
        raise ValueError('When tracking duplicate fragments, no reads were found in input BAM. Check input BAM: %s' % args.input_bam)

    # Finally, write out duplicate counts and insert sizes
    if args.duplicate_read_counts:
        with open(args.duplicate_read_counts, 'w') as output_file:
            output_file.write('\t'.join(['cell', 'total', 'total_deduplicated']) + '\n')

            for cell in cell_duplicate_counts:
                deduplicated_total_fragments, total_fragments = cell_duplicate_counts[cell]

                # To keep with convention of rest of pipeline, convert fragments to read counts (paired end)
                deduplicated_total_reads = deduplicated_total_fragments * 2
                total_reads = total_fragments * 2
                output_file.write('\t'.join([cell, str(total_reads), str(deduplicated_total_reads)]) + '\n')

    if args.insert_sizes:
        with open(args.insert_sizes, 'w') as output_file:
            output_file.write('\t'.join(['insert_size', 'count']) + '\n')

            for insert_size,count in enumerate(insert_sizes):
                output_file.write('\t'.join([str(insert_size), str(count)]) + '\n')
