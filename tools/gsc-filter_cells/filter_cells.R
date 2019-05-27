# First step of the signature-based workflow
# Remove low quality cells below a user-fixed cutoff of 
# percentiles or raw values of number of genes detected or
# total aligned reads

# Example of command (that generates output files) :
# Rscript filter_cells.R -f <input> --sep "/t" --absolute_genes 1700 --absolute_counts 90000 --pdfplot <plotfile.pdf>  --output <filtered_cells.tsv> --output_metada <filtered_metadata.tsv>

# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)
library(ggplot2)

# Arguments
option_list = list(
  make_option(c("-f", "--file"), default=NA, type='character',
              help="Input file that contains values to filter"),
  make_option("--sep", default="/t", type='character',
              help="File column separator [default : '%default' ]"),
  make_option("--percentile_genes", default=0, type='integer',
              help="nth Percentile of the number of genes detected by a cell distribution (If --percentile TRUE) [default : '%default' ]"),
  make_option("--percentile_counts", default=0, type='integer',
              help="nth Percentile of the total counts per cell distribution (If --percentile TRUE) [default : '%default' ]"),
  make_option("--absolute_genes", default=0, type='integer',
              help="Remove cells that didn't express at least this number of genes (Used if --percentile FALSE) [default : '%default' ]"),
  make_option("--absolute_counts", default=0, type='integer',
              help="Number of transcript threshold for cell filtering (Used if --percentile FALSE) [default : '%default' ]"),
  make_option("--pdfplot", type = 'character',
              help="Path to pdf file of the plots"),
  make_option("--output", type = 'character',
              help="Path to tsv file of filtered cell data"),
  make_option("--output_metada", type = 'character',
              help="Path to tsv file of filtered cell metadata")
)
opt = parse_args(OptionParser(option_list=option_list), args = commandArgs(trailingOnly = TRUE))

# check consistency of filtering options
if ((opt$percentile_counts > 0) & (opt$absolute_counts > 0)) {
  opt$percentile_counts = 100 } # since input parameters are not consistent (one or either method, not both), no filtering
if ((opt$percentile_counts == 0) & (opt$absolute_counts == 0)) {
  opt$percentile_counts = 100 } # since input parameters are not consistent (one or either method, not both), no filtering
if ((opt$percentile_genes > 0) & (opt$absolute_genes > 0)) {
  opt$percentile_genes = 100 } # since input parameters are not consistent (one or either method, not both), no filtering
if ((opt$percentile_genes == 0) & (opt$absolute_genes == 0)) {
  opt$percentile_genes = 100 } # since input parameters are not consistent (one or either method, not both), no filtering

# Import datasets
data.counts <- read.table(
  opt$file,
  header = TRUE,
  stringsAsFactors = F,
  sep = opt$sep,
  check.names = FALSE,
  row.names = 1
)

QC_metrics <-
  data.frame(cell_id = colnames(data.counts),
             nGenes = colSums(data.counts != 0),    # nGenes : Number of detected genes for each cell
             total_counts = colSums(data.counts),  # total_counts : Total counts per cell
             stringsAsFactors = F)

plot_hist <- function(mydata, variable, title, cutoff){
  mybinwidth = round(max(mydata[, variable]) * 5 / 100)
  mylabel = paste0("cutoff= ", cutoff)
  hist_plot <- qplot(
    mydata[, variable],
    main = title,
    xlab = variable,
    geom="histogram",
    binwidth = mybinwidth,
    col = I("white")) +
    geom_vline(xintercept = cutoff) +
    annotate(geom="text",
             x=cutoff + mybinwidth, y=1,
             label=mylabel, color="white")
    plot(hist_plot)
}

# returns the highest value such as the sum of the ordered values including this highest value
# is lower (below) than the percentile threshold (n)
percentile_cutoff <- function(n, qcmetrics, variable, plot_title, ...){
  p = n / 100
  percentile_threshold = quantile(qcmetrics[, variable], p)[[1]]
  plot_hist(qcmetrics,
            variable,
            plot_title,
            percentile_threshold)
  return(percentile_threshold)
}

pdf(file = opt$pdfplot)

# Determine thresholds based on percentile

if (opt$percentile_counts > 0) {
  counts_threshold <- percentile_cutoff(
    opt$percentile_counts,
    QC_metrics,
    "total_counts",
    "Histogram of Aligned read counts per cell"
  )} else {
  counts_threshold <- opt$absolute_counts
  plot_hist(QC_metrics,
            variable = "total_counts",
            title = "Histogram of Total counts per cell",
            cutoff = counts_threshold)
}

if (opt$percentile_genes > 0) {
 
  genes_threshold <- percentile_cutoff(
    opt$percentile_genes,
    QC_metrics,
    "nGenes",
    "Histogram of Number of detected genes per cell"
  )} else {
  genes_threshold <- opt$absolute_genes
  plot_hist(QC_metrics,
            variable = "nGenes",
            title = "Histogram of Number of detected genes per cell",         
            cutoff = genes_threshold)
}

# Filter out rows below thresholds (genes and read counts)
QC_metrics$filtered <- (QC_metrics$nGenes < genes_threshold) | (QC_metrics$total_counts < counts_threshold)

## Plot the results

# Determine title from the parameter logics
if(opt$percentile_counts > 0){
    part_one = paste0("Cells with aligned reads counts below the ",
                      opt$percentile_counts,
                      "th percentile of aligned read counts")} else {
    part_one = paste0("Cells with aligned read counts below ",
                      opt$absolute_counts)
}

if(opt$percentile_genes > 0){
    part_two = paste0("cells with number of detected genes below the ",
                      opt$percentile_genes,
                      "th percentile of detected gene counts.")} else {
    part_two = paste0("cells with number of detected genes below ",
                      opt$absolute_genes)
}

# plot with ggplot2
ggplot(QC_metrics, aes(nGenes, total_counts, colour = filtered)) +
     geom_point() + scale_y_log10() +
     scale_colour_discrete(name  = "",
                           breaks= c(FALSE, TRUE),
                           labels= c(paste0("Not filtered (", table(QC_metrics$filtered)[1], " cells)"), 
                                     paste0("Filtered (", table(QC_metrics$filtered)[2], " cells)"))) +
     xlab("Detected genes per cell") + ylab("Aligned reads per cell (log10 scale)") +
     geom_vline(xintercept = genes_threshold) + geom_hline(yintercept = counts_threshold) +
     ggtitle( paste0(part_one, " and\n", part_two, "\nwere filtered out")) +
     theme(plot.title = element_text(size = 8, face = "bold"))

dev.off()

# Retrieve identifier of kept cells
kept.cells <- QC_metrics$cell_id[!QC_metrics$filtered]

data.counts <- data.frame(Genes=rownames(data.counts[,kept.cells]), data.counts[,kept.cells])

# Save filtered cells 
write.table(
  data.counts,
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

# Add QC metrics of filtered cells to a metadata file
metadata <- QC_metrics[kept.cells,]

# Save the metadata (QC metrics) file
write.table(
  metadata,
  opt$output_metada,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)
