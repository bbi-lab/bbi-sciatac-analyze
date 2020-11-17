library(monocle3)
library(ggplot2)
library(readr)
library(argparse)
library(stringr)
library(Matrix)

#monocle_version = packageVersion("monocle3")
#print('Checking Monocle version...')
#if (seurat_version < 3) {
#    stop("Monocle 3 or greater is required due totype breaking changes in their API.")
#}
library(argparse)
library(Matrix)

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
source(paste0(dirname(thisFile()), '/', 'r_helper_functions/dim_reduction.R'))

get_aux_files = function(mtx_file) {
	base_name = str_replace(mtx_file, '[.]gz$', '')
	base_name = str_replace(base_name, '[.]mtx$', '')

	features_file = paste0(base_name, '.rows.txt')
	cells_file = paste0(base_name, '.columns.txt')

	return(list(features=features_file, cells=cells_file))
}

get_cell_metadata <- function(cells_file) {
	cell_metadata<-read.table(cells_file, header=FALSE, sep='\t', stringsAsFactors=FALSE)
    names(cell_metadata)[1]<-'V1'
    row.names(cell_metadata)<-cell_metadata$V1
	return(cell_metadata)
}

get_feature_metadata <- function(features_file) {
	feature_metadata<-data.frame(read.table(dim_file$features,header=FALSE, sep='\t', stringsAsFactors=FALSE))
	feature_metadata<-data.frame(feature_metadata)
	names(feature_metadata)<-c('gene_id', 'gene_short_name', 'biotype')
	row.names(feature_metadata)<-feature_metadata$gene_id
	return(feature_metadata )
}

plot_umap_pdf <- function(in_cds, umap_plot_file) {
	umap_plot_file <- str_replace(umap_plot_file, '[.]pdf', '')
    umap_plot_file <- paste0(umap_plot_file, '.pdf')
	plot_cells(cds, reduction_method='UMAP', show_trajectory_graph=FALSE)
    ggsave(umap_plot_file)
}

parser <- ArgumentParser(description='Script to perform do TF-IDF-based dim reduction on sci-ATAC-seq data.')
parser$add_argument('peak_matrix_file', help='MTX file to load.')
parser$add_argument('promoter_matrix_file', help='MTX file to load with promoter counts.')
parser$add_argument('--pca_coords', required=TRUE, help='PCA coords.')
parser$add_argument('--umap_coords', required=TRUE, help='UMAP coords.')
parser$add_argument('--tsne_coords', required=TRUE, help='TSNE coords.')
parser$add_argument('--tfidf_matrix', required=TRUE, help='TFIDF matrix.')
parser$add_argument('--umap_plot', required=TRUE, help='UMAP plot PDF file')
parser$add_argument('--monocle3_cds', help='Monocle3 CDS')
parser$add_argument('--svd_dimensions', default=75, type='integer', help='Max dimensions to calculate in SVD step.')
# parser$add_argument('--include_svd_1', action='store_true', help='By default the first dim is dropped from SVD as others have recommended in literature. Overall not make huge diff. Set this flag to use all dims.')
parser$add_argument('--sites_per_cell_threshold', default=100, type='integer', help='The min number of sites nonzero in a cell allowed. Cells not meeting this are filtered.')
parser$add_argument('--remove_top_ntile', type='double', help='Remove the top specified fraction of cells by total count. 0.975 would remove the top 2.5% of cells.')
parser$add_argument('--fast_tsne_path', help='Path to FAST TSNE executable.')
parser$add_argument('--regress_depth', action='store_true', help='Regress depth from PC space using limma.')
parser$add_argument('--exclude_chromosomes', nargs='+', help='List of chromosomes to exclude from dimensionality reduction. Sex chromosomes could be a good example of this.')
args <- parser$parse_args()

# if (args$include_svd_1) {
#     svd_start_dim <- 1    
# } else {
#     svd_start_dim <- 2
# }

# Load data
print('Loading data...')
promoter_matrix <- load_mtx_file(args$promoter_matrix_file)
peak_matrix <- load_mtx_file(args$peak_matrix_file)

# Binarize matrix
peak_matrix@x[peak_matrix@x > 0] <- 1

# Filter any chromosomes as needed
if (! is.null(args$exclude_chromosomes)) {
    invalid_peaks <- grepl(paste(paste0(args$exclude_chromosomes, '_'), collapse="|"), 
                        rownames(peak_matrix))
    print('filtered requested chromosomes')
    print(table(invalid_peaks))
    peak_matrix <- peak_matrix[!invalid_peaks, ]
}

valid_cells <- intersect(colnames(promoter_matrix), colnames(peak_matrix))

if (length(valid_cells) == 0) {
    warning('No cells overlap between peak matrix and promoter_matrix provided.')
    quit( save='no', status=0 )
}

peak_matrix <- peak_matrix[, valid_cells]
promoter_matrix <- promoter_matrix[, valid_cells]

# Filter peak matrix
print('Filtering input data...')
peak_matrix <- filter_features(peak_matrix, cells=ncol(peak_matrix) * 0.01)
peak_matrix <- filter_cells(peak_matrix, args$sites_per_cell_threshold)

## Remove upper outliers if requested
total_counts <- Matrix::colSums(peak_matrix)

if (!is.null(args$remove_top_ntile)) {
    upper_threshold <- quantile(total_counts, 1 - args$remove_top_ntile)
    peak_matrix <- peak_matrix[, total_counts < upper_threshold]
    total_counts <- total_counts[total_counts < upper_threshold]
}
filtered_cells <- colnames(peak_matrix)

if (length(filtered_cells) == 0) {
    warning('No peak matrix cells remain after filtering.')
    quit( save='no', status=0 )
}

promoter_matrix <- promoter_matrix[, colnames(peak_matrix)]

#print('Computing TFIDF...')
#tf_idf_counts <- tfidf(peak_matrix)
##
#print('Computing SVD...')
#cell_embeddings <- do_svd(tf_idf_counts, dims=args$svd_dimensions)
#
#if (args$regress_depth) {
#    regress_from_pca(cell_embeddings, log10(total_counts))
#}

print('Creating Monocle3 CDS...')
dim_file <- get_aux_files(args$promoter_matrix_file)

cell_metadata <- get_cell_metadata(dim_file$cells)
cell_metadata<-cell_metadata[cell_metadata$V1 %in% filtered_cells,,drop=FALSE]

feature_metadata <- get_feature_metadata(dim_file$features)

cds <- new_cell_data_set(expression_data=promoter_matrix,
                         cell_metadata=cell_metadata,
                         gene_metadata=feature_metadata)

print('Pre-processing CDS...')
cds <- detect_genes(cds)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, method="LSI")

print('Reducing dimensions...')
cds <- reduce_dimension(cds, reduction_method='tSNE', preprocess_method='LSI')
cds <- reduce_dimension(cds, reduction_method='UMAP', preprocess_method='LSI')

pca_coords <- reducedDims(cds)$LSI
tsne_coords <- reducedDims(cds)$tSNE
umap_coords <- reducedDims(cds)$UMAP

print('Clustering cells...')
cds <- cluster_cells(cds, reduction_method='UMAP', resolution=1.0e-5)

print('Writing result...')
#write_mtx_file(tf_idf_counts, args$tfidf_matrix)

readr::write_delim(data.frame(umap_coords), file=args$umap_coords, delim='\t')
readr::write_delim(data.frame(tsne_coords), file=args$tsne_coords, delim='\t')
readr::write_delim(data.frame(pca_coords), file=args$pca_coords, delim='\t')

plot_umap_pdf(in_cds=cds, args$umap_plot)

if (!is.null(args$monocle3_cds)) {
    saveRDS(cds, args$monocle3_cds)

}
