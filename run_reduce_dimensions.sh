#!/bin/bash

SAMPLE="sample1"
#
# Variables
#
# SCRIPTS_DIR="/net/gs/vol1/home/bge/src/andrew_hill/sciatac_pipeline/src"
SCRIPTS_DIR="/net/gs/vol1/home/bge/src/sciatac_pipeline/sciatac_pipeline-master/src"
MATRIX_DIR="/net/gs/vol1/home/bge/cole/sci-atac-seq/tests/gregs_three_level/analysis_run1/make_matrices"
REDUCED_DIR="/net/gs/vol1/home/bge/cole/sci-atac-seq/tests/gregs_three_level/analysis_run1/reduce_dimension"
FAST_TSNE_PATH="$SCRIPTS_DIR/FIt-SNE/bin/fast_tsne"

#
# Files
#
PEAK_MATRIX="$MATRIX_DIR/sample1.peak_matrix.mtx.gz"
PROMOTER_MATRIX="$MATRIX_DIR/sample1.promoter_matrix.mtx.gz"
SVD_COORDS="$REDUCED_DIR/$SAMPLE.svd_coords.txt"
UMAP_COORDS="$REDUCED_DIR/$SAMPLE.umap_coords.txt"
TSNE_COORDS="$REDUCED_DIR/$SAMPLE.tsne_coords.txt"
TFIDF_MATRIX="$REDUCED_DIR/$SAMPLE.tfidf_matrix.mtx"

SVD_DIMENSIONS="75"
SITES_PER_CELL_THRESHOLD="100"
REMOVE_TOP_NTILE="0.025"

# SEURAT_OBJECTS="$REDUCED_DIR/$SAMPLE.seurat_object.rds"

module load zlib/1.2.6 pigz/latest


Rscript $SCRIPTS_DIR/reduce_dimensions.R \
        $PEAK_MATRIX \
        $PROMOTER_MATRIX \
        --svd_coords $SVD_COORDS \
        --umap_coords $UMAP_COORDS \
        --tsne_coords $TSNE_COORDS \
        --tfidf_matrix $TFIDF_MATRIX \
        --svd_dimensions $SVD_DIMENSIONS \
        --sites_per_cell_threshold $SITES_PER_CELL_THRESHOLD \
        --remove_top_ntile $REMOVE_TOP_NTILE \
        --fast_tsne_path $FAST_TSNE_PATH

