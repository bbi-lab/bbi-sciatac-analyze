library(Matrix)
library(cisTopic)
library(stringr)
library(argparse)
library(varhandle)
library(S4Vectors)
library(GenomicRanges)
library(rGREAT)
library(RcisTarget)
library(Rtsne)
library(monocle)
library(plyr)
library(dplyr)

thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
source(paste0(dirname(thisFile()), '/', 'r_helper_functions/io_functions.R'))


parser = argparse::ArgumentParser(description='Script to generate a cisTopicModel for a given subset of the mouse data.')
parser$add_argument('matrix', help='Binary matrix from mouse sci-ATAC-seq dataset.')
parser$add_argument('model', help='RDS output file containing model.')
parser$add_argument('--metadata', help='Table of cell metadata. Must contain a "cell" column.')
parser$add_argument('--field', help='Field in mouse metadata to subset on.')
parser$add_argument('--value', help='Value within metadata column specified by --field to subset.')
parser$add_argument('--site_percent_min', default=0.005, type='double', help='Min percent of cells where site can be detected.')
parser$add_argument('--min_cell_nonzero', default=100, type='integer', help='Min number of non-zero entries that a cell can have within the matrix. Cells not meeting are filtered.')
parser$add_argument('--topics', nargs="+", required=TRUE, help='List of topic numbers to run.')
parser$add_argument('--no_great', action='store_true', help='Set flag to turn off GREAT enrichments.')
parser$add_argument('--promoter_matrix', action='store_true', help='Set flag to indicate that the matrix is a promoter matrix and genomic ranges should not be expected as names for rows.')
args = parser$parse_args()

args$topics = as.numeric(args$topics)

print('Loading data...')
if (str_detect(args$matrix, '[.]rds$')) {
	atac_matrix = readRDS(args$matrix)
} else if (str_detect(args$matrix, '[.]mtx')) {
	atac_matrix = load_mtx_file(args$matrix)
} else {
	stop('Matrix format not supported. Only .rds or .mtx files are allowed.')
}

# Binarize (if not already)
atac_matrix@x[atac_matrix@x > 0] = 1

if (!is.null(args$metadata) & !is.null(args$value) & !is.null(args$column)) {
	metadata = read.delim(args$metadata)

	if (varhandle::check.numeric(args$value)) {
		args$value = as.numeric(args$value)
	}

	if (! args$field %in% colnames(metadata)) {
		stop('--field argument specified is not a column in metadata...')
	}
	cells = metadata[metadata[, args$field] == args$value, ]$cell
} else {
	cells = colnames(atac_matrix)
}

if (! length(cells) > 1) {
	stop('1 or fewer cells found with specified --value in --field column. Try a larger group of cells / check you did not mess anything up.')
}

print('Subsetting data...')
valid_cells = cells[Matrix::colSums(atac_matrix[, cells] > 0) > args$min_cell_nonzero]
atac_matrix.subset = atac_matrix[, valid_cells]
site_totals = Matrix::rowSums(atac_matrix.subset)
atac_matrix.subset = atac_matrix.subset[site_totals > (length(cells) * args$site_percent_min) & site_totals < (length(cells) * args$site_percent_max), ]

if (!args$promoter_matrix) {
	coords = str_split_fixed(rownames(atac_matrix.subset), '_', 3)
	new_coords = paste0(coords[, 1], ':', coords[, 2], '-', coords[, 3])
	rownames(atac_matrix.subset) = new_coords
} else {
	rownames(atac_matrix.subset) = paste0('chr1', ':', 1:nrow(atac_matrix.subset), '-', 2:(nrow(atac_matrix.subset) + 1))
}

# Go through the topic modeling vignette (WARNING LONG RUNTIMES LIKELY)
print('Training models...')
cisTopicObject <- createcisTopicObject(atac_matrix.subset, project.name='mouse_atac', keepCountsMatrix = FALSE)
cisTopicObject <- runModels(cisTopicObject, topic=args$topics, seed=2018, nCores=1, burnin = 250, iterations = 500)
print('Done with models.')

cisTopicObject = selectModel(cisTopicObject)
cisTopicObject = getRegionsScores(cisTopicObject, method='Zscore', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=FALSE)

cisTopicObject = cisTopic::runtSNE(cisTopicObject, perplexity=30)

saveRDS(cisTopicObject, args$model)

if (!args$no_great) {
	cisTopicObject <- GREAT(cisTopicObject, genome='hg38', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
	saveRDS(cisTopicObject, args$model)
}
