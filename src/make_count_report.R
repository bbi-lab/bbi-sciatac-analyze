options(stringsAsFactors=FALSE)
library(argparse)
library(readr)
library(dplyr)

parser = argparse::ArgumentParser(description='Script to merge count tables produced by the pipeline into one single table.')
parser$add_argument('duplicate_counts', help='Duplicate count table.')
parser$add_argument('peak_counts', help='Table of counts from peak regions.')
parser$add_argument('tss_counts', help='Table of counts from tss regions.')
parser$add_argument('count_report', help='Output file with all counts merged.')
args = parser$parse_args()

# duplicate counts file format
# cell    total   total_deduplicated
# P2-E12_P2-E12_F09-rowF09-colF06 19840   12748
duplicate_counts = readr::read_delim(args$duplicate_counts, '\t')

# peak counts file format
# cell    count
# P1-E01_P1-E01_A01-rowE01-colA05 1846
peak_counts = readr::read_delim(args$peak_counts, '\t')

# tss counts file format
# cell    count
# P1-E01_P1-E01_A01-rowE01-colA05 996
tss_counts = readr::read_delim(args$tss_counts, '\t')

colnames(peak_counts) = c('cell', 'total_deduplicated_peaks')
colnames(tss_counts) = c('cell', 'total_deduplicated_tss')

count_report = duplicate_counts %>%
    dplyr::left_join(peak_counts, by='cell') %>%
    dplyr::left_join(tss_counts, by='cell')

count_report$total_deduplicated_peaks[is.na(count_report$total_deduplicated_peaks)] = 0
count_report$total_deduplicated_tss[is.na(count_report$total_deduplicated_tss)] = 0

readr::write_delim(count_report, path=args$count_report, delim='\t')
