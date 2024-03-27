#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(monocle3)
  library(argparse)
})


createCDSWithHash <- function(in_cds_path, merged_hash_filename) {
    # Load the cds object of peaks by cells
    timeseriesPilotFull <- readRDS(in_cds_path)
    # Load the hash dataframe
    mergedHashInfo <- readRDS(merged_hash_filename)
    rownames(mergedHashInfo) <- mergedHashInfo$CellID

    # Incorporate the hash information as columns
    mergedColData <- merge(colData(timeseriesPilotFull), mergedHashInfo, by = 0, all = FALSE)
    rownames(mergedColData) <- mergedColData$cell
    cdsObjWithHash <- timeseriesPilotFull[,rownames(colData(timeseriesPilotFull)) %in% rownames(mergedColData)]
    # The rownames match between the two objects now
    # rownames(mergedColData) == rownames(colData(mergedSetOfColumns))
    colData(cdsObjWithHash) <- mergedColData
    return(cdsObjWithHash)
}

# createCDSWithHash <- function(dirToCDS,pathToHash,prefixName){
#     # Load the cds object of peaks by cells
#     fileCDS = paste0(prefixName,'-monocle3_cds.rds')
#     pathToCDS <- paste0(dirToCDS,fileCDS)
#     timeseriesPilotFull <- readRDS(pathToCDS)
#     # Load the hash dataframe
#     mergedHashInfo <- readRDS(paste0(pathToHash,'/mergedHashInfo.RDS'))
#     rownames(mergedHashInfo) <- mergedHashInfo$CellID
# 
#     # Incorporate the hash information as columns
#     mergedColData <- merge(colData(timeseriesPilotFull), mergedHashInfo, by = 0, all = FALSE)
#     rownames(mergedColData) <- mergedColData$cell
#     cdsObjWithHash <- timeseriesPilotFull[,rownames(colData(timeseriesPilotFull)) %in% rownames(mergedColData)]
#     # The rownames match between the two objects now
#     # rownames(mergedColData) == rownames(colData(mergedSetOfColumns))
#     colData(cdsObjWithHash) <- mergedColData
#     return(cdsObjWithHash)
# }

parser <- argparse::ArgumentParser(description = "Script to add ATAC-seq hash read information to a CDS colData.")
parser$add_argument("--in_cds_name", help = "An input cell_data_set filename.")
parser$add_argument("--in_hash_info_rds_name", help = "An input hash information RDS filename.")
parser$add_argument("--out_cds_name", help = "An output CDS filename.")
args <- parser$parse_args()

out_cds <- createCDSWithHash(in_cds_path=args$in_cds_name, merged_hash_filename=args$in_hash_info_rds_name)
saveRDS(out_cds, file=args$out_cds_name)

