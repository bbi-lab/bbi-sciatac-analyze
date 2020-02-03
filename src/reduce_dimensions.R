library(Seurat)
library(readr)
library(argparse)
library(stringr)

seurat_version = packageVersion("Seurat")
print('Checking seurat version...')
if (seurat_version < 3) {
    stop("Seurat 3 or greater is required due breaking changes in their API.")
}
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

get_embeddings_df = function(seurat_obj, embedding) {
    embeddings = Seurat::Embeddings(object=seurat_obj, reduction=embedding)
    cells = as.vector(rownames(embeddings))

    embeddings = as.data.frame(embeddings)
    dims = colnames(embeddings)
    embeddings$cell = cells

    embeddings = embeddings[, c('cell', dims)]
    return(embeddings)
}

parser = ArgumentParser(description='Script to perform do TF-IDF-based dim reduction on sci-ATAC-seq data.')
parser$add_argument('peak_matrix_file', help='MTX file to load.')
parser$add_argument('promoter_matrix_file', help='MTX file to load with promoter counts.')
parser$add_argument('--svd_coords', required=TRUE, help='TSNE coords.')
parser$add_argument('--umap_coords', required=TRUE, help='UMAP coords.')
parser$add_argument('--tsne_coords', required=TRUE, help='TSNE coords.')
parser$add_argument('--tfidf_matrix', required=TRUE, help='TSNE coords.')
parser$add_argument('--seurat_object', help='TSNE coords.')
parser$add_argument('--svd_dimensions', default=75, type='integer', help='Max dimensions to calculate in SVD step.')
parser$add_argument('--include_svd_1', action='store_true', help='By default the first dim is dropped from SVD as others have recommended in literature. Overall not make huge diff. Set this flag to use all dims.')
parser$add_argument('--sites_per_cell_threshold', default=100, type='integer', help='The min number of sites nonzero in a cell allowed. Cells not meeting this are filtered.')
parser$add_argument('--remove_top_ntile', type='double', help='Remove the top specified fraction of cells by total count. 0.975 would remove the top 2.5% of cells.')
parser$add_argument('--fast_tsne_path', help='Path to FAST TSNE executable.')
parser$add_argument('--regress_depth', action='store_true', help='Regress depth from PC space using limma.')
parser$add_argument('--exclude_chromosomes', nargs='+', help='List of chromosomes to exclude from dimensionality reduction. Sex chromosomes could be a good example of this.')
args = parser$parse_args()

if (args$include_svd_1) {
    svd_start_dim = 1    
} else {
    svd_start_dim = 2
}

# Load data
print('Loading data...')
promoter_matrix = load_mtx_file(args$promoter_matrix_file)
peak_matrix = load_mtx_file(args$peak_matrix_file)

# Binarize matrix
peak_matrix@x[peak_matrix@x > 0] = 1

# Filter any chromsomes as needed
if (! is.null(args$exclude_chromosomes)) {
    invalid_peaks = grepl(paste(paste0(args$exclude_chromosomes, '_'), collapse="|"), 
                        rownames(peak_matrix))
    print('filtered requested chromosomes')
    print(table(invalid_peaks))
    peak_matrix = peak_matrix[!invalid_peaks, ]
}

valid_cells = intersect(colnames(promoter_matrix), colnames(peak_matrix))

if (length(valid_cells) == 0) {
    stop('No cells overlap between peak matrix and promoter_matrix provided.')
}

peak_matrix = peak_matrix[, valid_cells]
promoter_matrix = promoter_matrix[, valid_cells]

# Filter peak matrix
print('Filtering input data...')
peak_matrix = filter_features(peak_matrix, cells=ncol(peak_matrix) * 0.01)
peak_matrix = filter_cells(peak_matrix, args$sites_per_cell_threshold)

## Remove upper outliers if requested
total_counts = Matrix::colSums(peak_matrix)

if (!is.null(args$remove_top_ntile)) {
    upper_threshold = quantile(total_counts, 1 - args$remove_top_ntile)
    peak_matrix = peak_matrix[, total_counts < upper_threshold]
    total_counts = total_counts[total_counts < upper_threshold]
}

promoter_matrix = promoter_matrix[, colnames(peak_matrix)]

print('Computing TFIDF...')
tf_idf_counts = tfidf(peak_matrix)

print('Computing SVD...')
cell_embeddings = do_svd(tf_idf_counts, dims=args$svd_dimensions)

if (args$regress_depth) {
    regress_from_pca(cell_embeddings, log10(total_counts))
}

print('Creating Seurat object...')
seurat_obj = Seurat::CreateSeuratObject(promoter_matrix)
seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', assay='RNA')

seurat_obj = Seurat::L2Dim(seurat_obj, reduction='pca')

print('Running UMAP...')
seurat_obj = Seurat::RunUMAP(seurat_obj, reduction = 'pca.l2', dims = svd_start_dim:args$svd_dimensions)

print('Computing tSNE...')
if (! is.null(args$fast_tsne_path)) {
    print('   using fast tSNE...')
    seurat_obj = Seurat::RunTSNE(seurat_obj,
                                reduction = 'pca.l2',
                                dims = svd_start_dim:args$svd_dimensions,
                                tsne.method = "FIt-SNE",
                                nthreads = 4,
                                reduction.name = "tsne",
                                reduction.key = "tSNE_",
                                fast_tsne_path = args$fast_tsne_path,
                                max_iter = 2000)
} else {
    print('   using Rtsne...')
    seurat_obj = Seurat::RunTSNE(seurat_obj,
                                 reduction = 'pca.l2',
                                 dims = svd_start_dim:args$svd_dimensions)
}

# Cluster cells at a few resolutions
seurat_obj = Seurat::FindNeighbors(seurat_obj, reduction='pca.l2', nn.eps=0.25, dims=svd_start_dim:args$svd_dimensions)
seurat_obj = Seurat::FindClusters(seurat_obj, reduction='pca.l2', n.start=20, resolution=0.1, dims=svd_start_dim:args$svd_dimensions)
seurat_obj = Seurat::FindClusters(seurat_obj, reduction='pca.l2', n.start=20, resolution=0.3, dims=svd_start_dim:args$svd_dimensions)
seurat_obj = Seurat::FindClusters(seurat_obj, reduction='pca.l2', n.start=20, resolution=0.6, dims=svd_start_dim:args$svd_dimensions)

# Write out results
write_mtx_file(tf_idf_counts, args$tfidf_matrix)

print('Saving results...')
readr::write_delim(get_embeddings_df(seurat_obj, 'pca'), path=args$svd_coords, delim='\t')
readr::write_delim(get_embeddings_df(seurat_obj, 'umap'), path=args$umap_coords, delim='\t')
readr::write_delim(get_embeddings_df(seurat_obj, 'tsne'), path=args$tsne_coords, delim='\t')

if (!is.null(args$seurat_object)) {
    saveRDS(seurat_obj, args$seurat_object)
}
