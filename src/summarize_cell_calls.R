library(rjson)
library(Matrix)
library(argparse)
library(readr)
library(dplyr)
library(stringr)
library(glue)
options(bitmapType='cairo')
options(stringsAsFactors = FALSE)

DOUBLET_PERCENTAGE_THRESHOLD = 0.9

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


#This code is directly adapted from the picard source code: http://picard-tools.sourcearchive.com/documentation/1.25-1/classnet_1_1sf_1_1picard_1_1sam_1_1DuplicationMetrics_a0a76ff106a01cb46966cbdf0d778843.html#a0a76ff106a01cb46966cbdf0d778843
fer = function(x,c,n){
  return(c/x - 1 + exp(-n/x))
}

picardcomp = function(c = 50,n = 100,m = 1, M = 10){
  if(c >= n || fer(m*c,c,n) <= 0){
    print("Invalid inputs!")
    break
  }
  while(fer(M*c,c,n) >= 0){
    M = M*10
  }
  
  for(i in 1:40){
    r = (m+M)/2
    u = fer(r*c,c,n)
    if(u == 0){
      break
    }else{
      if(u > 0){
        m = r
      }else{
        if(u < 0){
          M = r
        }
      }
    }
  }
  return(c*(m+M)/2)
}

bloom_collision <- function(n1, n2, n12){
  #from https://github.com/jbloomlab/multiplet_freq
  n <- n1 * n2 / n12
  mu1 <- -log((n - n1) / n)
  mu2 <- -log((n - n2) / n)
  mu <- mu1 + mu2
  return (1 - mu * exp(-mu) / (1 - exp(-mu)))
}

parser = argparse::ArgumentParser(description='Script to make summary plots of cell calling results and other basic summary stats.')
parser$add_argument('--stats_files', nargs='+', required=TRUE, help='List of all JSON files containing cell calling stats.')
parser$add_argument('--sample_names', nargs='+', required=TRUE, help='Sample names.')
parser$add_argument('--read_count_tables', nargs='+', required=TRUE, help='Original read count tables (not thresholded).')
parser$add_argument('--insert_size_tables', nargs='+', required=TRUE, help='Insert size distribution tables.')
parser$add_argument('--peak_groups', nargs='+', required=TRUE, help='Peak groups assigned to sample.')
parser$add_argument('--peak_files', nargs='+', required=TRUE, help='External peak files assigned to sample.')
parser$add_argument('--peak_call_files', nargs='+', required=TRUE, help='BED files with peak calls.')
parser$add_argument('--merged_peaks', required=TRUE, help='BED file with merged peak set.')
parser$add_argument('--per_base_tss_region_coverage_files', nargs='+', required=TRUE, help='Set of files with per base coverage for TSS regions.')
parser$add_argument('--combined_duplicate_report', nargs='+', required=TRUE, help='Combined non-mitochondrial and mitochondrial duplicate report files.')
parser$add_argument('--window_matrices', nargs='+', required=FALSE, help='Set of files with per base coverage for TSS regions.')
parser$add_argument('--barnyard', action='store_true', help='Set if sample is a barnyard sample.')
parser$add_argument('--plot', required=TRUE, help='Plot summarizing results.')
parser$add_argument('--output_stats', required=TRUE, help='Output file containing stats displayed in the plots.')
args = parser$parse_args()

total_merged_peaks = nrow(read.delim(args$merged_peaks, header=FALSE))

pdf(args$plot,height=11,width=8)
par(mfrow=c(3,2))

# TODO Using an lapply here to plot AND recover the output stats is pretty terrible...
# but I could not get the for loop way (storing in list as go) to work in a rush... woo R!
output_stats = lapply(1:length(args$stats_files), function(i) {
  sample_stats = rjson::fromJSON(file=args$stats_files[[i]])
#   sample_name = stringr::str_replace(basename(args$stats_files[[i]]), '.stats.json', '')
  sample_name = args$sample_names[[i]]
  peak_group = args$peak_groups[[i]]
  peak_file = args$peak_files[[i]]
  message(glue('Processing {sample_name}...'))

  # read count table file format
  # cell    total   total_deduplicated      total_deduplicated_peaks        total_deduplicated_tss
  # P2-E12_P2-E12_F09-rowF09-colF06 19840   12748   3942    1858
  message('-> loading counts...')
  sample_counts = readr::read_delim(args$read_count_tables[[i]], '\t')

  # combined read count table file format (includes mitochondrial read columns)
  # cell    total_nonmito   total_nonmito_deduplicated        total_mito      total_mito_deduplicated
  # P2-E12_P2-E12_F09-rowF09-colF06 19840   12748   404     246
  message('-> loading combined read counts...')
  combined_read_counts = readr::read_delim(args$combined_duplicate_report[[i]], '\t')

  # insert sizes table file format
  # insert_size     count
  # 0       0
  # 1       0
  message('-> loading insert sizes...')
  sample_insert_sizes = readr::read_delim(args$insert_size_tables[[i]], '\t')

  message('-> loading sample peaks...')
  sample_peak_counts <- tryCatch(
                                  {
                                    nrow(read.delim(args$peak_call_files[[i]], header=FALSE))
                                  },
                                  error = function(cond) {
                                    return(0)
                                  },
                                  warning = function(cond) {
                                    return(0)
                                  }
                                )

  # tss region coverage file format
  # position        total
  # -1000.5 1534
  message('-> loading TSS region coverage...')
  sample_tss_coverage = readr::read_delim(args$per_base_tss_region_coverage_files[[i]], '\t')

  if (args$barnyard) {
    message('-> loading window matrix to calculate doublet stats...')
    sample_window_matrix = load_mtx_file(args$window_matrices[[i]])
    human_counts = Matrix::colSums(sample_window_matrix[grepl('HUMAN_', rownames(sample_window_matrix)),])
    mouse_counts = Matrix::colSums(sample_window_matrix[grepl('MOUSE_', rownames(sample_window_matrix)),])

    if (sum(human_counts) == 0) {
      stop('No human counts found. Chromosomes must start with HUMAN_ in name (e.g. HUMAN_19) for human and MOUSE_ for mouse.')
    }
    if (sum(mouse_counts) == 0) {
      stop('No mouse counts found. Chromosomes must start with MOUSE_ in name (e.g. MOUSE_19) for mouse and HUMAN_ for human.')
    }

    barnyard_df = data.frame(cell=names(mouse_counts), human=human_counts, mouse=mouse_counts)

    # Cells with > 90% of one genome vs. other are singlets
    human_cell = with(barnyard_df, human/(human + mouse) > DOUBLET_PERCENTAGE_THRESHOLD)
    mouse_cell = with(barnyard_df, mouse/(human + mouse) > DOUBLET_PERCENTAGE_THRESHOLD)

    # Doublets get colored in red
    barnyard_df$color = 'black'
    barnyard_df$color[!(human_cell | mouse_cell)] = 'red'

    mousecounts = sum(mouse_cell)
    humancounts = sum(human_cell)
    totalcounts = nrow(barnyard_df)
    humanfrac = humancounts/(humancounts + mousecounts)
    mousefrac = mousecounts/(humancounts + mousecounts)
    collisioninflation = ((humanfrac^2) + (mousefrac^2))/(2*humanfrac*mousefrac)
    bloom_collision_rate = signif(bloom_collision(totalcounts-mousecounts,totalcounts-humancounts,totalcounts-(humancounts+mousecounts)),4)
  }

  # Read distributions
  currcellfloor = sample_stats$cell_threshold
  currsubcells = sample_counts$total_deduplicated >= currcellfloor
  
  fraction_hs = round(sum(sample_counts[currsubcells, 'total_deduplicated_peaks'])/sum(sample_counts[currsubcells,"total_deduplicated"]),3)
  fraction_tss = round(sum(sample_counts[currsubcells, 'total_deduplicated_tss'])/sum(sample_counts[currsubcells,"total_deduplicated"]),3)
  total_reads = sum(sample_counts$total)
  total_deduplicated_reads <- sum(sample_counts$total_deduplicated)
  fraction_reads_in_cells = sum(sample_counts[currsubcells, 'total']) / total_reads
  total_barcodes = nrow(sample_counts)
  number_of_cells = nrow(sample_counts[currsubcells, ])
  median_reads_per_cell = median(sample_counts$total_deduplicated[currsubcells])
  min_reads_per_cell = min(sample_counts$total_deduplicated[currsubcells])
  max_reads_per_cell = max(sample_counts$total_deduplicated[currsubcells])

  if(length(currsubcells) == 0){next}

  count_data = sample_counts$total_deduplicated + 1
  rank_data = rank(-count_data, ties.method='random')
  ymin = 1
  ymax = max(count_data) * 1.1
  xmin = 1
  xmax = max(rank_data)

  ## Show first 100, last 100 for noise and cells
  ## and then show 2000 (or less) evenly spaced points for each.
  ## This helps keep the plot sizes smaller.
  cells_at_border = 500
  cell_borders = c(seq(1, cells_at_border), seq(number_of_cells - cells_at_border, number_of_cells))
  cell_interval = max(floor(number_of_cells / 2000), 1)
  cell_points =  c(cell_borders, seq(1, number_of_cells, by=cell_interval))
  
  noise_border = c(seq(number_of_cells + 1, number_of_cells + cells_at_border), seq(total_barcodes - cells_at_border, total_barcodes))
  noise_interval = max(floor((total_barcodes - number_of_cells) / 2000), 1)
  noise_points = c(noise_border, seq(number_of_cells + 1, total_barcodes, by=noise_interval))

  cell_points = rank_data %in% cell_points
  noise_points = rank_data %in% noise_points
  
  cell_df = data.frame(x=rank_data[cell_points], y=count_data[cell_points])
  noise_df = data.frame(x=rank_data[noise_points], y=count_data[noise_points])

  ## Now make the plot
  # plot: knee plot
  plot(x=cell_df$x, y=cell_df$y, col="mediumseagreen", main=paste0(sample_name, ' log10(Reads Per Cell)'),
       ylab="Number of Reads", xlab='Barcode Rank', pch=16, log='xy', cex=1, ylim = c(ymin, ymax), xlim=c(xmin, xmax))
  points(x=noise_df$x, y=noise_df$y, col='#d3d3d3', log='xy', pch=16, cex=1, ylim = c(ymin, ymax), xlim=c(xmin, xmax))


  legend("bottomleft",c(paste0("Total Reads: ", total_reads),
                    paste0("\n Fraction Reads in Cells: ", round(fraction_reads_in_cells, 4)),
                    paste0("\n Total Barcodes: ",total_barcodes),
                    paste0("\n Number of Cells: ", number_of_cells),
                    paste0("\n Fraction HS (cells only): ", round(fraction_hs, 4)),
                    paste0("\n Median Reads/Cell: ", median_reads_per_cell),
                    paste0("\n Range of Reads/Cell: ", min_reads_per_cell," - ", max_reads_per_cell)),bty="n", cex=0.75, pt.cex = 1, text.font=2)
  grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")

  file_name <- paste0(sample_name, '-knee_plot.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  plot(x=cell_df$x, y=cell_df$y, col="mediumseagreen", main=paste0(sample_name, ' log10(Reads Per Cell)'),
       ylab="Number of Reads", xlab='Barcode Rank', pch=16, log='xy', cex=1, ylim = c(ymin, ymax), xlim=c(xmin, xmax))
  points(x=noise_df$x, y=noise_df$y, col='#d3d3d3', log='xy', pch=16, cex=1, ylim = c(ymin, ymax), xlim=c(xmin, xmax))
  legend("bottomleft",c(paste0("Total Reads: ", total_reads),
                    paste0("\n Fraction Reads in Cells: ", round(fraction_reads_in_cells, 4)),
                    paste0("\n Total Barcodes: ",total_barcodes),
                    paste0("\n Number of Cells: ", number_of_cells),
                    paste0("\n Fraction HS (cells only): ", round(fraction_hs, 4)),
                    paste0("\n Median Reads/Cell: ", median_reads_per_cell),
                    paste0("\n Range of Reads/Cell: ", min_reads_per_cell," - ", max_reads_per_cell)),bty="n", cex=0.75, pt.cex = 1, text.font=2)
  grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")  
  dev.off()

  # Most other plots use only called cells
  subsetmat.sample = sample_counts[currsubcells, ]

  # Duplicate reporting
  # plot: fragment distribution
  deduped = which(subsetmat.sample$total > subsetmat.sample$total_deduplicated)
  totalfrags = apply(subsetmat.sample[deduped,c("total", "total_deduplicated")],1,function(x){picardcomp(x[2],x[1])})

  median_duplication_rate = median(with(subsetmat.sample, (total - total_deduplicated) / total))
  median_fraction_molecules_observed = median(subsetmat.sample$total_deduplicated[deduped]/totalfrags)
  median_total_fragments = median(totalfrags)

  hist(log10(totalfrags), main=paste0(sample_name, " Estimated Fragments"), col="orange",lwd=2,pch=20,las=1, breaks=60)
  legend("topright",
    c(paste0("\n Median Frac Mol Obs: ", round(median_fraction_molecules_observed, 4)),
    paste0("\n Median Dup Rate: ", round(median_duplication_rate, 4)),
    paste0("\n Median Total Frags Est: ", round(median_total_fragments, 4))),
  bty="n", cex=0.75, pt.cex = 1, text.font=2)
  abline(v=log10(median_total_fragments),lwd=2,lty="dashed")

  file_name <- paste0(sample_name, '-estimated_frag_dist.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  hist(log10(totalfrags), main=paste0(sample_name, " Estimated Fragments"), col="orange",lwd=2,pch=20,las=1, breaks=60)
  legend("topright",
    c(paste0("\n Median Frac Mol Obs: ", round(median_fraction_molecules_observed, 4)),
    paste0("\n Median Dup Rate: ", round(median_duplication_rate, 4)),
    paste0("\n Median Total Frags Est: ", round(median_total_fragments, 4))),
  bty="n", cex=0.75, pt.cex = 1, text.font=2)
  abline(v=log10(median_total_fragments),lwd=2,lty="dashed")
  dev.off()
  
  # FRiP plot
  # plot: frip
  frip = subsetmat.sample$total_deduplicated_peaks/subsetmat.sample$total_deduplicated
  median_per_cell_frip = median(frip)

  hist(frip,
        xlab="Fraction of Reads Mapping to DHS", main=paste0(sample_name, ' Cells FRiP'), col="dodgerblue2",lwd=2,pch=20,las=1, breaks=60)
  abline(v=median_per_cell_frip,lwd=2,lty="dashed")
  legend("topright",
    c(paste0("\n Median FRiP: ", round(median_per_cell_frip, 4)),
    paste0("\n Sample Peaks: ", sample_peak_counts),
    paste0("\n Total Merged Peaks: ", total_merged_peaks),
    paste0("\n Peak group: ", peak_group),
    paste0("\n Peak file: ", peak_file)),
  bty="n", cex=0.75, pt.cex = 1, text.font=2)

  file_name <- paste0(sample_name, '-frip.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  hist(frip,
        xlab="Fraction of Reads Mapping to DHS", main=paste0(sample_name, ' Cells FRiP'), col="dodgerblue2",lwd=2,pch=20,las=1, breaks=60)
  abline(v=median_per_cell_frip,lwd=2,lty="dashed")
  legend("topright",
    c(paste0("\n Median FRiP: ", round(median_per_cell_frip, 4)),
    paste0("\n Sample Peaks: ", sample_peak_counts),
    paste0("\n Total Merged Peaks: ", total_merged_peaks),
    paste0("\n Peak group: ", peak_group),
    paste0("\n Peak file: ", peak_file)),
  bty="n", cex=0.75, pt.cex = 1, text.font=2)
  dev.off()
  
  # TSS plot
  # plot: frit
  frit = subsetmat.sample$total_deduplicated_tss/subsetmat.sample$total_deduplicated
  median_per_cell_frit= median(frit)
  hist(frit,
        xlab="Fraction of Reads Mapping to TSS Regions", main=paste0(sample_name, ' Cells FRiT'), col="red",lwd=2,pch=20,las=1, breaks=60)
  abline(v=median_per_cell_frit,lwd=2,lty="dashed")
  legend("topright",
    c(paste0("\n Median FRiT: ", round(median_per_cell_frit, 4))),
    bty="n", cex=0.75, pt.cex = 1, text.font=2)

  file_name <- paste0(sample_name, '-frit.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  hist(frit,
        xlab="Fraction of Reads Mapping to TSS Regions", main=paste0(sample_name, ' Cells FRiT'), col="red",lwd=2,pch=20,las=1, breaks=60)
  abline(v=median_per_cell_frit,lwd=2,lty="dashed")
  legend("topright",
    c(paste0("\n Median FRiT: ", round(median_per_cell_frit, 4))),
    bty="n", cex=0.75, pt.cex = 1, text.font=2)
  dev.off()

  # Insert size distribution
  # plot: insert sizes
  plot(x=sample_insert_sizes$insert_size, y=sample_insert_sizes$count, type='l', col='red', xlab='Insert Size (bp)', ylab='Count (reads)', main=paste0(sample_name, ' Insert Sizes'))
  grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")
  
  file_name <- paste0(sample_name, '-insert_size_dist.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  plot(x=sample_insert_sizes$insert_size, y=sample_insert_sizes$count, type='l', col='red', xlab='Insert Size (bp)', ylab='Count (reads)', main=paste0(sample_name, ' Insert Sizes'))
  grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")
  dev.off()
  
  # TSS enrichment  
  ## compute the min at edge of window as mean of left and right for average of 10bp on each side
  left_count = mean(sample_tss_coverage$total[sample_tss_coverage$position < min(sample_tss_coverage$position) + 10])
  right_count = mean(sample_tss_coverage$total[sample_tss_coverage$position > max(sample_tss_coverage$position) - 10])
  min_count = mean(c(left_count, right_count))

  ## compute enrichment over min and then define tss enrichment as max in the 200bp preceeding the center of the window
  sample_tss_coverage$enrichment = sample_tss_coverage$total / min_count
  tss_enrichment = max(sample_tss_coverage$enrichment[sample_tss_coverage$position > -200 & sample_tss_coverage$position < 0]) 

  # For barnyard sub out TSS enrichment plot for barnyard plot
  if (args$barnyard) {
  	# plot: barnyard
    plot(barnyard_df$human,barnyard_df$mouse,pch=20,xlab="Human reads",ylab="Mouse reads", col=barnyard_df$color)
    abline(a=0, b=1-DOUBLET_PERCENTAGE_THRESHOLD,lwd=2,lty="dashed", col='lightgrey')
    abline(a=0, b=1/(1-DOUBLET_PERCENTAGE_THRESHOLD),lwd=2,lty="dashed", col='lightgrey')

    legend("topright",c(paste0("Total Cells: ",totalcounts),
                        paste0("Human Cells: ",humancounts),
                        paste0("Mouse Cells: ",mousecounts),
                        paste0("Observed Collision Rate: ",signif(1-(humancounts + mousecounts)/totalcounts,4)),
                        paste0("Calculated Collision Rate: ",signif(2*collisioninflation*(1-(humancounts + mousecounts)/totalcounts),4)),
                        paste0("Bloom Collision Rate: ",bloom_collision_rate)))

    file_name <- paste0(sample_name, '.png')
    png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
    plot(barnyard_df$human,barnyard_df$mouse,pch=20,xlab="Human reads",ylab="Mouse reads", col=barnyard_df$color)
    abline(a=0, b=1-DOUBLET_PERCENTAGE_THRESHOLD,lwd=2,lty="dashed", col='lightgrey')
    abline(a=0, b=1/(1-DOUBLET_PERCENTAGE_THRESHOLD),lwd=2,lty="dashed", col='lightgrey')

    legend("topright",c(paste0("Total Cells: ",totalcounts),
                        paste0("Human Cells: ",humancounts),
                        paste0("Mouse Cells: ",mousecounts),
                        paste0("Observed Collision Rate: ",signif(1-(humancounts + mousecounts)/totalcounts,4)),
                        paste0("Calculated Collision Rate: ",signif(2*collisioninflation*(1-(humancounts + mousecounts)/totalcounts),4)),
                        paste0("Bloom Collision Rate: ",bloom_collision_rate)))
    dev.off()
  } else {
  	# plot: tss_enrichment
    plot(sample_tss_coverage$position, sample_tss_coverage$enrichment, type='l', col='black', xlab='Position relative to TSS (bp)', ylab='Fold Enrichment', main=paste0(sample_name, ' TSS enrichment'))
    grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")
    legend("topright",
        c(paste0("\n TSS enrichment: ", round(tss_enrichment, 4))),
    bty="n", cex=0.75, pt.cex = 1, text.font=2)

    file_name <- paste0(sample_name, '-tss_enrichment.png')
    png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
    plot(sample_tss_coverage$position, sample_tss_coverage$enrichment, type='l', col='black', xlab='Position relative to TSS (bp)', ylab='Fold Enrichment', main=paste0(sample_name, ' TSS enrichment'))
    grid(nx = 10, ny = 5, col = "lightgray", lty = "dotted")
    legend("topright",
        c(paste0("\n TSS enrichment: ", round(tss_enrichment, 4))),
    bty="n", cex=0.75, pt.cex = 1, text.font=2)
    dev.off()
  }

  # UMI per cell histogram
  # plot: read distribution (umi per cell)
  file_name <- paste0(sample_name, '-umi_per_cell.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  x <- sample_counts$total_deduplicated
  x <- x[x>10]
  xrange <- range(log10(x))
  breaks <- 20 * (xrange[2] - xrange[1])
  hist(log10(x),xaxt='n',axes=TRUE,breaks=breaks,plot=TRUE,main='Histogram of UMI per Cell',xlab='UMI per cell (log scale)')
  tick_location <- axTicks(1)
  tick_value <- 10^tick_location
  axis(1, at=tick_location, labels=tick_value)
  dev.off()
 
  # Mitochondrial read fraction histogram.
  subsetmitomat.sample = combined_read_counts[combined_read_counts$total_nonmito_deduplicated >= currcellfloor, ]
  fraction_mitochondrial_reads <- subsetmitomat.sample$total_mito_deduplicated/(subsetmitomat.sample$total_nonmito_deduplicated+subsetmitomat.sample$total_mito_deduplicated)
  hist(fraction_mitochondrial_reads, main=paste0(sample_name," Fraction of Mitochondrial Reads"), col="darkorchid1", lwd=2,pch=20,las=1, breaks=60)
  if(sum(subsetmitomat.sample$total_mito_deduplicated) == 0) {
    plot_lims <- par('usr')
    text((plot_lims[2]-plot_lims[1])*0.5+plot_lims[1], (plot_lims[4]-plot_lims[3])*0.5+plot_lims[3], 'No mitochondrial reads.')
  }
  file_name <- paste0(sample_name, '-fraction_mitochondrial_reads.png')
  png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
  hist(fraction_mitochondrial_reads, main=paste0(sample_name," Fraction of Mitochondrial Reads"), col="darkorchid1", lwd=2,pch=20,las=1, breaks=60)
  if(sum(subsetmitomat.sample$total_mito_deduplicated) == 0) {
    plot_lims <- par('usr')
    text((plot_lims[2]-plot_lims[1])*0.5+plot_lims[1], (plot_lims[4]-plot_lims[3])*0.5+plot_lims[3], 'No mitochondrial reads.')
  }
  dev.off()
 
  message('-> sample done.')
  if (args$barnyard) {
    stats_df <- data.frame('sample'=sample_name,
                           'cell_threshold'=currcellfloor,
                           'fraction_hs'=fraction_hs,
                           'fraction_tss'=fraction_tss,
                           'median_per_cell_frip'=median_per_cell_frip,
                           'median_per_cell_frit'=median_per_cell_frit,
                           'tss_enrichment'=tss_enrichment,
                           'sample_peaks_called'=sample_peak_counts,
                           'total_merged_peaks'=total_merged_peaks,
                           'total_reads'=total_reads,
                           'fraction_reads_in_cells'=fraction_reads_in_cells,
                           'total_barcodes'=total_barcodes,
                           'number_of_cells'=number_of_cells,
                           'median_reads_per_cell'=median_reads_per_cell,
                           'min_reads_per_cell'=min_reads_per_cell,
                           'max_reads_per_cell'=max_reads_per_cell,
                           'median_duplication_rate'=median_duplication_rate,
                           'median_fraction_molecules_observed'=median_fraction_molecules_observed,
                           'median_total_fragments'=median_total_fragments,
                           'total_deduplicated_reads'=total_deduplicated_reads,
                           'bloom_collision_rate'=bloom_collision_rate)
  } else {
    stats_df <- data.frame('sample'=sample_name,
                           'cell_threshold'=currcellfloor,
                           'fraction_hs'=fraction_hs,
                           'fraction_tss'=fraction_tss,
                           'median_per_cell_frip'=median_per_cell_frip,
                           'median_per_cell_frit'=median_per_cell_frit,
                           'tss_enrichment'=tss_enrichment,
                           'sample_peaks_called'=sample_peak_counts,
                           'total_merged_peaks'=total_merged_peaks,
                           'total_reads'=total_reads,
                           'fraction_reads_in_cells'=fraction_reads_in_cells,
                           'total_barcodes'=total_barcodes,
                           'number_of_cells'=number_of_cells,
                           'median_reads_per_cell'=median_reads_per_cell,
                           'min_reads_per_cell'=min_reads_per_cell,
                           'max_reads_per_cell'=max_reads_per_cell,
                           'median_duplication_rate'=median_duplication_rate,
                           'median_fraction_molecules_observed'=median_fraction_molecules_observed,
                           'median_total_fragments'=median_total_fragments,
                           'total_deduplicated_reads'=total_deduplicated_reads)
  }
  # Return any key stats for output file
  return(stats_df)
})
dev.off()

output_df = bind_rows(output_stats)
readr::write_delim(output_df, path=args$output_stats, delim='\t')
