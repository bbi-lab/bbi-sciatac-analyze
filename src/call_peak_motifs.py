# This code is adapated from the ANNOTATE_PEAKS stage of 10X scATAC pipeline
# Most of the functions here are taken with very minor modification from
# cellranger-atac. The interface, etc. obviously changed to fit in with our pipeline.
from __future__ import print_function
import numpy as np
import tempfile
import pyfaidx
import MOODS
import MOODS.scan
import MOODS.tools
import MOODS.parsers
import argparse

from Bio import motifs
from collections import namedtuple, Counter, OrderedDict

JASPAR_MOTIF_BASE_ORDER = ['A', 'C', 'G', 'T']
ALLOWED_FORMATS = ['bed', 'motif', 'position', 'count', 'binary', 'binary-bed']


class Motifs:

    def __init__(self, fasta, motifs_input, bg=None):
        self.all_motifs = []
        with open(motifs_input, "r") as infile:
            self.all_motifs = list(motifs.parse(infile, "jaspar"))

        # for large sequence header, only keep the text before the first space
        self.genome_seq = pyfaidx.Fasta(fasta, one_based_attributes=False, key_function=lambda x: x.split()[0])
        self.bg = bg

    def scan_motif_from_bed(self, peaks_iter, tf_genes=None, out_file=None, out_format='bed', use_genome_bg=True,
                            pseudocount=1.0, pvalue=1e-5, window=7):
        """
        :param peaks_iter: peaks iterator, yielding a tuple (chrom, start, end, strand)
        :param tf_genes: a list of gene symbols for motif search. If None, scan all motifs in the collection
        :param out_file: file name of output. Default is no output file
        :param out_format: output format can be bed, position, count or binary
        :param use_genome_bg: whether to generate background model using the genome seq.
                              If False, sequences from the input bed file will be used for background
        :param pseudocount: MOODS scan pseudocounts for the calculation of log-odds scores from matrices
        :param pvalue: motif hits cutoff
        :param window: the scanning window size of the lookahead filtration algorithm
        :return:
        """

        assert out_format in ALLOWED_FORMATS
        if use_genome_bg:
            self.bg = self.get_reference_bg(self.genome_seq)

        motif = self.get_motif_of(tf_genes)

        # Each TF gets a matrix for the + and for the - strand, and a corresponding threshold
        (matrices, thresholds) = self._prepare_moods_settings(motif, self.bg, pseudocount, pvalue)

        scanner = MOODS.scan.Scanner(window)  # parameter is the window size
        scanner.set_motifs(matrices, self.bg, thresholds)

        if out_file is not None:
            out = open(out_file, 'w')

        maps = []
        for peak_idx, peak in enumerate(peaks_iter):

            bed_coord = {'chr': peak[0],
                         'start': int(peak[1]),
                         'stop': int(peak[2]),
                         'name': peak_idx,
                         'strand': peak[3]}
            # original: pyfasta.Fasta access
#            seq = self.genome_seq.sequence(bed_coord, one_based=False)
            # pyfaidx.Fasta access
            if(peak[3] in (-1, '-1', '-')):
                seq = self.genome_seq[peak[0]][int(peak[1]):int(peak[2])].reverse.complement
            else:
                seq = self.genome_seq[peak[0]][int(peak[1]):int(peak[2])]

            # seq is of unicode format, need to convert to str
            results = scanner.scan(str(seq))
            parsed = self._parse_scan_results(results, motif, bed_coord, out_format)

            if out_file is not None:
                out.writelines(['\t'.join(map(str, item)) + '\n' for item in parsed])
            else:
                maps.extend(parsed)

        if out_file is not None:
            out.close()
        else:
            return maps

    def get_motif_of(self, tf_genes):

        # tf_genes is None, return all motifs in the collection
        if tf_genes is None:
            return self.all_motifs

        selected = []
        for motif in self.all_motifs:
            motif_name = motif.name.split('_')[0]
            if motif_name in tf_genes:
                selected.append(motif)
        return selected

    @staticmethod
    def get_reference_bg(fa):
        b = bytearray()
        for chrom in fa.keys():
            b.extend(bytes(fa[chrom]))

        print("Computing the background model...")
        bg = MOODS.tools.bg_from_sequence_dna(str(b), 1)
        print("The background model is {}\n".format(bg))

        return bg

    def _prepare_moods_settings(self, jaspar_motifs, bg, pseduocount, pvalue=1e-5):
        """Find hits of list of jaspar_motifs in pyfasta object fasta, using the background distribution bg and
        pseudocount, significant to the give pvalue
        """
        matrices = [self._jaspar_to_moods_matrix(j, bg, pseduocount) for j in jaspar_motifs]
        matrices = matrices + [MOODS.tools.reverse_complement(m) for m in matrices]
        thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]
        return matrices, thresholds

    @staticmethod
    def _jaspar_to_moods_matrix(jaspar_motif, bg, pseudocount):
        """Convert JASPAR motif into a MOODS log_odds matrix, using the give background distribution & pseudocounts
        - JASPAR 2010 matrix_only format::

          >MA0001.1 AGL3
          A  [ 0  3 79 40 66 48 65 11 65  0 ]
          C  [94 75  4  3  1  2  5  2  3  3 ]
          G  [ 1  0  3  4  1  0  5  3 28 88 ]
          T  [ 2 19 11 50 29 47 22 81  1  6 ]

        - JASPAR 2010-2014 PFMs format::

          >MA0001.1 AGL3
          0       3       79      40      66      48      65      11      65      0
          94      75      4       3       1       2       5       2       3       3
          1       0       3       4       1       0       5       3       28      88
          2       19      11      50      29      47      22      81      1       6

        """
        with tempfile.NamedTemporaryFile() as fn:
            f = open(fn.name, "w")
            for base in 'ACGT':
                line = " ".join(str(x) for x in jaspar_motif.pwm[base])
                f.write(line + "\n")

            f.close()
            return MOODS.parsers.pfm_to_log_odds(fn.name, bg, pseudocount)

    @staticmethod
    def _parse_scan_results(moods_scan_res, motifs, bed_coord, out_format="bed"):
        """ Parse results of MOODS.scan.scan_dna and return/write
            The default input contains one pair of a single motif: forward and reverse strand

            out_format: if "bed" each hit will be the motif position in bed format, with the start and end as peak coordinates
                        if "motif" same as "bed" except the start and end columns are the motif location
                        if "position" each hit will be a list of positions within a peak region (relative to the start position)
                        if "count" each hit will be an integer of the number of occurrences. NO OUTPUT if count == 0
                        if "binary" each hit will be 1 (any occurrence).  NO OUTPUT if count == 0
        """
        assert out_format in ALLOWED_FORMATS

        all_hits = []
        for (motif_idx, hits) in enumerate(moods_scan_res):

            motif = motifs[motif_idx % len(motifs)]
            strand = "-" if motif_idx >= len(motifs) else "+"

            if len(hits) > 0:

                if out_format == "binary":
                    # for binary format we only care about whether len(hits)>0
                    record = [motif_idx % len(motifs), bed_coord['name'], 1 if len(hits) > 0 else 0]
                    all_hits.append(record)
                    continue

                elif out_format == "count":
                    # for count format we just need len(hits)
                    record = [motif_idx % len(motifs), bed_coord['name'], len(hits)]
                    all_hits.append(record)
                    continue

                elif out_format == "binary-bed":
                    record = [bed_coord['chr'], bed_coord['start'], bed_coord['stop'], motif.name]
                    all_hits.append(record)
                    continue

                else:
                    # for bed or position format, we ouput each individual hit
                    # output bed file of [chr start end motifname score strand pos pos+motiflen]

                    for h in hits:
                        if out_format == "bed":
                            motif_start = bed_coord['start'] + int(h.pos)
                            motif_end = bed_coord['start'] + int(h.pos) + motif.length
                            score = round(h.score, 4)
                            record = [bed_coord['chr'], bed_coord['start'], bed_coord['stop'], motif.name, score, strand, motif_start, motif_end]

                        elif out_format == "motif":
                            motif_start = bed_coord['start'] + int(h.pos)
                            motif_end = bed_coord['start'] + int(h.pos) + motif.length
                            score = round(h.score, 4)
                            record = [bed_coord['chr'], motif_start, motif_end, motif.name, score, strand]

                        else:
                            strand = 1 if strand == "-" else -1
                            record = [motif_idx % len(motifs), bed_coord['name'], int(h.pos), round(h.score, 4), strand]

                        all_hits.append(record)

        return all_hits

def combine_csv(input_csvs, output_csv, header_lines=1):
    """Combine a list of CSV files specified in input_csvs into
    a single output csv in output_csv. It is assumed that all
    CSVs have the same header structure. The number of header
    lines is specified in header_lines
    """
    if not input_csvs:
        output_csv = None
        return
    with open(output_csv, "w") as out:
        for i, icsv in enumerate(input_csvs):
            with open(icsv, "r") as infile:
                header = []
                for h in xrange(header_lines):
                    header.append(infile.next())
                if i == 0:
                    for hline in header:
                        out.write(hline)
                for line in infile:
                    out.write(line)

                    
Peak = namedtuple('Peak', ['chrom', 'start', 'end', 'strand'])
def peak_reader(peaks, select=[]):
    '''Streaming peak reader with selection mode. Accepts ordered indices'''
    curr = 0
    with open(peaks, 'r') as peak_file:
        for line, row in enumerate(peak_file):
            peak = row.strip("\n").split("\t")
            if curr < len(select) and line == select[curr]:
                curr += 1
                continue
            yield Peak(peak[0], int(peak[1]), int(peak[2]), '.')

def get_peak_GC_counts(peak, genome_fa, counts=True):
    '''Get GC% in base seq in a peak and (optionally) nucleotide counts'''

    seq = str(genome_fa[peak.chrom][peak.start:peak.end]).upper()
    base_counter = {}
    base_counter['A'] = seq.count('A') + seq.count('a')
    base_counter['G'] = seq.count('G') + seq.count('g')
    base_counter['C'] = seq.count('C') + seq.count('c')
    base_counter['T'] = seq.count('T') + seq.count('t')
    peakGC = (base_counter['G'] + base_counter['C']) / sum(base_counter.values())
    if counts:
        return peakGC, base_counter
    else:
        return peakGC


def get_GCbinned_peaks_and_bg(peaks, genome_fa, GCbins, pseudocount=1.0):
    '''Calculate the G-C and A-T % within the peaks for each bin'''

    def findbin(peakGC, GCbins):
        '''Expected sorted GCbins'''
        for gc in GCbins:
            if peakGC < gc[1] and peakGC >= gc[0]:
                return gc
        if abs(peakGC - GCbins[-1][1]) < 1e-5:
            return GCbins[-1]

    GCdict = OrderedDict()
    base_counter = {}
    for gc in GCbins:
        GCdict[gc] = {}
        GCdict[gc]['counter'] = Counter({base: 0 for base in JASPAR_MOTIF_BASE_ORDER})
        GCdict[gc]['peaks'] = []

    for num, peak in enumerate(peaks):
        peakGC, base_counter = get_peak_GC_counts(peak, genome_fa)
        gc = findbin(peakGC, GCbins)
        GCdict[gc]['peaks'].append(num)
        GCdict[gc]['counter'] += Counter(base_counter)

    for gc in GCbins:
        GCdict[gc]['total_bases'] = sum(GCdict[gc]['counter'].values())

        freq = {}
        for base in JASPAR_MOTIF_BASE_ORDER:
            freq[base] = np.divide(GCdict[gc]['counter'][base] + pseudocount,
                                   GCdict[gc]['total_bases'] + 4 * pseudocount, dtype=float)

        # the freq is computed using only the "+" strand. To get the freq of both strand, average A&T and G&C
        bg = {}
        bg['A'] = (freq['A'] + freq['T']) / 2
        bg['C'] = (freq['G'] + freq['C']) / 2
        bg['G'] = bg['C']
        bg['T'] = bg['A']
        GCdict[gc]['counter'] = [bg[base] for base in JASPAR_MOTIF_BASE_ORDER]

    return GCdict

def load_gc_bins(fasta, peaks):
    genome_fa = pyfaidx.Fasta(fasta, one_based_attributes=True, key_function=lambda x: x.split()[0])

    # get peak-GC distribution
    GCdist = [get_peak_GC_counts(peak, genome_fa, counts=False) for peak in peak_reader(peaks)]

    # compute base background from peaks in bins
    GCbounds = np.percentile(GCdist, np.linspace(0, 100, 26, endpoint=True), interpolation='lower')
    GCbins = sorted(list(set(zip(GCbounds, GCbounds[1:]))))  # uniqify
    peaks = peak_reader(peaks)
    GCdict = get_GCbinned_peaks_and_bg(peaks, genome_fa, GCbins)
    return GCdict


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to generate motif calls within a set of peaks binned by GC content.')
    parser.add_argument('fasta', help='Fasta file to use as reference.')
    parser.add_argument('peaks', help='Peaks in BED format to use as reference.')
    parser.add_argument('motifs', help='Motif database (JASPAR format) to use as motif reference.')
    parser.add_argument('output_file', help='Output BED file with motif calls.')
    parser.add_argument('--pwm_threshold', default=1e-7, type=float, help='Threshold to consider a PWM hit a hit.')
    parser.add_argument('--gc_bin', required=True, type=int, help='Which GC bin (1 to 25) to limit processing to.')
    args = parser.parse_args()

    gc_bins = load_gc_bins(args.fasta, args.peaks)
    
    gc_bin = list(gc_bins.keys())[args.gc_bin]
    low, high = gc_bin
    print('processing bin from GC content %s to %s...' % (low, high))

    motifs = Motifs(args.fasta, args.motifs, bg=gc_bins[gc_bin]['counter'])
    peaks_iter = peak_reader(args.peaks, select=gc_bins[gc_bin]['peaks'])
        
    motifs.scan_motif_from_bed(peaks_iter,
                               out_file=args.output_file,
                               out_format="binary-bed",
                               use_genome_bg=False,
                               pvalue=args.pwm_threshold)
