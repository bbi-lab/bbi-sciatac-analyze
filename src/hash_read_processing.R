## From Andrew Mullen in slack DM channel 20240111 with minor modifications (bge) 20240223

## Andrew Mullen
# 11/15/23
# The goal of this script is to plot the sciATAC-seq data and the hashing information

# This script has a variety of utility R functions that can be used to determine QC metrics from a hashRDS directory generated from the ATAC-seq pipeline

# it will return five figures
# 1) The total hash captured by each cell
# 2) The total of each hash well captured 
# 3) Top to Second Ratios w/o Normalization
# 4) Top to Second Ratios w/ Normalization for total hash measured
# 5) A barchart of cells that pass hashing 

# It also will return a dataframe of relevant hashing infromation

library(ggplot2)
library(tidyr)

# Need to generate a hashInfo.rds object that includes all of the relevant information from hashing


hashInfoRDSGen <- function(sampleName, called_cells_file, hash_umis_per_cell_file, output_plot_file="HistogramHashesCounter.pdf") {
    # Load the cellIDs that are getting called as cells during the pipeline
    cellsCalled = read.table(called_cells_file, header =TRUE, sep = "\t")
    # Get the hash reads
    hashesCalled = read.table(hash_umis_per_cell_file, header =FALSE, sep = "\t")
    mergeDFHashing <- hashesCalled[hashesCalled$V2 %in% cellsCalled$cell,]
    # Create a histogram of these hashes captured
    hashDF <- data.frame(hashCounter = mergeDFHashing$V3)
    rownames(hashDF) <- mergeDFHashing$V2
    histoGenForCellsCaptured(hashDF, output_plot_file)
    # Save the hashDF object
    return(hashDF)
}


# Make a histogram of the number of hash molecules captured per cell
# Function to make a histogram of the number of hash molecules captured per cell
# Inputs are hashDF - a dataframe of the cells that were called and the has counts associated with that cell
# output_plot_file - name of output PDF file.
histoGenForCellsCaptured <- function(hashDF, output_plot_file="HistogramHashesCounter.pdf") {
    plot <- ggplot(hashDF, aes(x = hashCounter)) +
    geom_histogram(binwidth =0.05,fill = "blue", color = "black") +
    scale_x_log10() +  # Set x-axis to log scale
    labs(x = "Log10 Number of Hash Molecules", y = "Number of Cells") +
    theme_minimal()
    # Save the plot as a PDF
    ggsave(output_plot_file, plot, width = 8, height = 6, units = "in")
}


# Count the number of hash molecules per hash to see if any biases so far
barchartHashCaptures <- function(hash_table_file="hashTable.out", output_plot_file="TotalHashesCaptured.pdf") {
    hashTable = read.table(hash_table_file, header =FALSE, sep = "\t")
    wide_tableHash <- spread(hashTable, key = V3, value = V5, fill = 0)
    rownames(wide_tableHash) <- wide_tableHash$V2
    wide_tableHashCorrected <- wide_tableHash[,4:99]    
    col_sums <- colSums(wide_tableHashCorrected)
    totalCapturedPerHash <- data.frame(col_sums)
    totalCapturedPerHash$HashName <- rownames(totalCapturedPerHash)
    totalHashes <- ggplot(totalCapturedPerHash, aes(x = HashName, y = col_sums)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        labs(title = "Hashes Captured By Well", x = "Hash Well", y = "Number of Hash") +
        theme_minimal()+
        theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))
        ggsave(output_plot_file, totalHashes, width = 12, height = 6, units = "in")
    return(wide_tableHashCorrected)
}


# Sort of the hashes to be organized A1 - H12
# Plot top to second ratios
find_largest_values <- function(row,colnames) {
    sorted_values <- sort(row, decreasing = TRUE)
    largest <- sorted_values[1]
    second_largest <- sorted_values[2]
    largest_index <- which(row == largest)
    hashRatio <- as.numeric(largest)/as.numeric(second_largest)
    largest_colname <- colnames[which(row == largest)[1]]
    return(list(largest, second_largest,hashRatio,largest_colname))
}


# This is taking a while to run, I think I should reduce the size of the matrix to only include the values that are in called cells
topToSecondHistogram <- function(hashMatrix, hashDF, output_plot_file="TopToSecondHashRatioLog10.pdf") {
    hashMatrixCellsOnly <- hashMatrix[rownames(hashMatrix) %in% rownames(hashDF),]
    hashTopToSecondTable <- apply(hashMatrixCellsOnly, 1, find_largest_values,colnames = colnames(hashMatrixCellsOnly))
    transposeHashTable <- do.call(rbind, hashTopToSecondTable)
    transposeHashTable <- data.frame(transposeHashTable)
    # Replace the NA values with values of 1000
    colnames(transposeHashTable) <- c("TopHash", "SecondHash", "ToptoSecondRatio","TopHashID")
    transposeHashTable$ToptoSecondRatio <- replace(unlist(transposeHashTable$ToptoSecondRatio), is.infinite(unlist(transposeHashTable$ToptoSecondRatio)), 10000)
    transposeHashTable$Log10TopToSecondRatio <- log10(transposeHashTable$ToptoSecondRatio)
    plot <- ggplot(transposeHashTable, aes(x = Log10TopToSecondRatio)) +
        geom_histogram(bins = 100,fill = "blue", color = "black") +  # Set x-axis to log scale
        labs(x = "Top To Second Ratio (Log 10)", y = "Count") +
        geom_vline(xintercept = log10(2), color = "red") +
        coord_cartesian(xlim = c(0, 2))      +
        theme_minimal()
    ggsave(output_plot_file, plot, width = 8, height = 6, units = "in")        
    return(transposeHashTable)
}


# Now do the same plot, but normalize based on the total hash captured
NormalizeHash <- function(hashMatrix, hashDF, output_plot_file="TopToSecondHashRatioNormLog10.pdf") {
    sumOfHashCaptured <- colSums(hashMatrix)
    hashMatNorm <- sweep(hashMatrix, MARGIN = 2, STATS = sumOfHashCaptured, FUN = "/")
    hashMatrixNormCellsOnly <- hashMatNorm[rownames(hashMatNorm) %in% rownames(hashDF),]
    hashTopToSecondTableNorm <- apply(hashMatrixNormCellsOnly, 1, find_largest_values,colnames = colnames(hashMatrixNormCellsOnly))
    transposeHashTableNorm <- do.call(rbind, hashTopToSecondTableNorm)
    transposeHashTableNorm <- data.frame(transposeHashTableNorm)    
    colnames(transposeHashTableNorm) <- c("TopHashNorm", "SecondHashNorm", "ToptoSecondRatioNorm","TopHashIDNorm")
    transposeHashTableNorm$ToptoSecondRatioNorm <- replace(unlist(transposeHashTableNorm$ToptoSecondRatioNorm), is.infinite(unlist(transposeHashTableNorm$ToptoSecondRatioNorm)), 10000)
    transposeHashTableNorm$ToptoSecondLog10RatioNorm <- log10(transposeHashTableNorm$ToptoSecondRatioNorm)
    plot <- ggplot(transposeHashTableNorm, aes(x = ToptoSecondLog10RatioNorm)) +
        geom_histogram(bins = 100,fill = "blue", color = "black") +  # Set x-axis to log scale
        labs(x = "Top To Second Ratio (Log 10)", y = "Count") +
        geom_vline(xintercept = log10(2), color = "red") +
        coord_cartesian(xlim = c(0, 2))        +
        theme_minimal()
    ggsave(output_plot_file, plot, width = 8, height = 6, units = "in")      
    return(transposeHashTableNorm)
}


# Make a funtion to plot how many cells are passing on hashQC
BarchartPassingHash <- function(hashTableTopToSecond, hashTableTopToSecondNorm, output_plot_file="mergeHashPassingHash.pdf", output_rds_file="mergedHashInfo.RDS") {
    # Form one dataframe with both
    hashTableTopToSecond$CellID <- rownames(hashTableTopToSecond)
    hashTableTopToSecondNorm$CellID <- rownames(hashTableTopToSecondNorm)
    mergeHashDF <- merge(hashTableTopToSecond,hashTableTopToSecondNorm,by="CellID",all=TRUE)
    passingHash <- (mergeHashDF$TopHash > 5) & (unlist(mergeHashDF$TopHashID) == unlist(mergeHashDF$TopHashIDNorm)) & (mergeHashDF$ToptoSecondRatioNorm > 2) & mergeHashDF$ToptoSecondRatio > 2
    passingHashNoMin <- (unlist(mergeHashDF$TopHashID) == unlist(mergeHashDF$TopHashIDNorm)) & (mergeHashDF$ToptoSecondRatioNorm > 2) & mergeHashDF$ToptoSecondRatio > 2
    mergeHashDF$PassingHash <- passingHash
    mergeHashDF$passingHashNoMin <- passingHashNoMin
    subsetDFForPassing <- mergeHashDF[,c(12,13)]
    dataFrameForPlotting = data.frame(passing = c(sum(mergeHashDF$PassingHash)/length(mergeHashDF$PassingHash),sum(mergeHashDF$passingHashNoMin)/length(mergeHashDF$PassingHash)))
    dataFrameForPlotting$failing = 1-dataFrameForPlotting$passing
    dataFrameForPlotting$condition = c('Min Of 5 Hash','No Min')
    meltedDF <- reshape2::melt(dataFrameForPlotting)
    p1 <- ggplot(meltedDF, aes(x = condition, y = value,fill=variable)) +
    geom_bar(stat = "identity",position = "dodge") +
    labs(title = "Percent of Cells Passing By Hash", x = "Categories", y = "Percentage") +
    theme_minimal()
    ggsave(output_plot_file,  p1, width = 8, height = 6, units = "in")
    percentPassing <- sum(passingHash)/dim(mergeHashDF)[1]
    saveRDS(mergeHashDF, output_rds_file)
    return(0)
}


# Make some QC plots that I like and an R .RDS file called mergedHashInfo.RDS which has the necessary info
createCDSWithHash <- function(input_cds_file="monocle3_cds.rds", input_hash_file="mergedHashInfo.RDS", prefixName, output_cds_file="mergedHashInfo.RDS") {
    # Load the cds object of peaks by cells
    timeseriesPilotFull <- readRDS(input_cds_file)
    # Load the hash dataframe
    mergedHashInfo <- readRDS(input_hash_file)
    rownames(mergedHashInfo) <- mergedHashInfo$CellID

    # Incorporate the hash information as columns
    mergedColData <- merge(colData(timeseriesPilotFull), mergedHashInfo, by = 0, all = FALSE)
    rownames(mergedColData) <- mergedColData$cell
    cdsObjWithHash <- timeseriesPilotFull[,rownames(colData(timeseriesPilotFull)) %in% rownames(mergedColData)]
    # The rownames match between the two objects now
    # rownames(mergedColData) == rownames(colData(mergedSetOfColumns))
    colData(cdsObjWithHash) <- mergedColData
    saveRDS(cdsObjWithHash)
    return(0)
}


######################################################################################
# main
######################################################################################

parser <- argparse::ArgumentParser(description='Script to plot the sciATAC-seq data and the hashing information.')
# Input arguments.
parser$add_argument('--sample_name', required=TRUE, help='Sample name.')
parser$add_argument('--called_cells_file', required=TRUE, help='File of called cells.')
parser$add_argument('--hash_table_file', required=TRUE, help='File of hash reads found.')
parser$add_argument('--hash_umis_per_cell_file', required=TRUE, help='File of UMIs per cell.')
# Output arguments.
parser$add_argument('--hashes_counter_histogram', required=TRUE, help='Output plot of hash counter.')
parser$add_argument('--total_hashes_captured_histogram', required=TRUE, help='Output plot of total hashes captured.')
parser$add_argument('--top_to_second_hash_ratio_histogram', required=TRUE, help='Output plot of top to second hash ratio.')
parser$add_argument('--normalized_top_to_second_hash_ratio_histogram',required=TRUE, help='Output plot of normalized top to second hash ratio.')
parser$add_argument('--merged_hash_passing_histogram', required=TRUE, help='Output plot of merged hash passing.')
parser$add_argument('--merged_hash_info', required=TRUE, help='Output file of merged hash info.')
args <- parser$parse_args()

# hash_df  <- hashInfoRDSGen(sampleName=args$sample_name, called_cells_file=args$called_cells_file,  hash_umis_per_cell_file=args$hash_umis_per_cell_file, output_plot_file="HistogramHashesCounter.pdf")
# hash_matrix  <- barchartHashCaptures(args$hash_table_file, output_plot_file="TotalHashesCaptured.pdf")
# hash_table_top_to_second  <- topToSecondHistogram(hash_matrix, hash_df, output_plot_file="TopToSecondHashRatioLog10.pdf")
# hash_table_top_to_second_norm  <- NormalizeHash(hash_matrix, hash_df, output_plot_file="TopToSecondHashRatioNormLog10.pdf")
# merged_hash_info  <- BarchartPassingHash(hash_table_top_to_second, hash_table_top_to_second_norm, output_plot_file="mergeHashPassingHash.pdf", output_rds_file="mergedHashInfo.RDS")

hash_df  <- hashInfoRDSGen(sampleName=args$sample_name, called_cells_file=args$called_cells_file,  hash_umis_per_cell_file=args$hash_umis_per_cell_file, output_plot_file=args$hashes_counter_histogram)
hash_matrix  <- barchartHashCaptures(args$hash_table_file, output_plot_file=args$total_hashes_captured_histogram)
hash_table_top_to_second  <- topToSecondHistogram(hash_matrix, hash_df, output_plot_file=args$top_to_second_hash_ratio_histogram)
hash_table_top_to_second_norm  <- NormalizeHash(hash_matrix, hash_df, output_plot_file=args$normalized_top_to_second_hash_ratio_histogram)
merged_hash_info  <- BarchartPassingHash(hash_table_top_to_second, hash_table_top_to_second_norm, output_plot_file=args$merged_hash_passing_histogram, output_rds_file=args$merged_hash_info)

