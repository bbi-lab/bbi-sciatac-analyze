library(dplyr)
library(readr)


parser = argparse::ArgumentParser(description='Script to aggregate results from per base coverage command.')
parser$add_argument('input_file', help='Output from command in PerBaseTSSRegionCoverage.')
parser$add_argument('output_file', help='Aggregated coverage across the window that accounts for strand.')
args = parser$parse_args()

# per base tss region coverage input format
# 1       3465586 3467587 ENSMUSG00000089699      1       +       54      1
# 1       3465586 3467587 ENSMUSG00000089699      1       +       55      1
col_names = c('chrom', 'start', 'end', 'gene', 'score', 'strand', 'position', 'count')
sample_tss_coverage = readr::read_delim(args$input_file, delim='\t', col_types='cddcdcdd', col_names=col_names)
sample_tss_coverage$position = sample_tss_coverage$position - 1 - max(sample_tss_coverage$position) / 2 # center region about zero

sample_tss_coverage$position[sample_tss_coverage$strand == '-'] = -sample_tss_coverage$position[sample_tss_coverage$strand == '-'] # reverse negative strand

per_position_aggregate = sample_tss_coverage %>%
    group_by(position) %>%
    summarize(total=sum(count))

readr::write_delim(per_position_aggregate, path=args$output_file, delim='\t')
