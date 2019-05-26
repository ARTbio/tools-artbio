#####################
#   Filter cells    #
#####################

# First step of the signature-based workflow
# Remove low quality cells bellow a user-fixed cutoff of 
# percentiles or raw values of number of genes detected or
# total aligned reads

# Example of command (that generates output files) :
# Rscript filter_cells.R -f <input> --sep "/t" --cutoff_genes 1700 --cutoff_counts 90000 --pdfplot <plotfile.pdf>  --output <filtered_cells.tsv> --output_metada <filtered_metadata.tsv>

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
  make_option("--cutoff_genes", default=2, type='integer',
              help="Remove cells that didn't express at least this number of genes (Used if --percentile FALSE) [default : '%default' ]"),
  make_option("--cutoff_counts", default=2, type='integer',
              help="Number of transcript threshold for cell filtering (Used if --percentile FALSE) [default : '%default' ]"),
  make_option("--pdfplot", type = 'character',
              help="Path to pdf file of the plots"),
  make_option("--output", type = 'character',
              help="Path to tsv file of filtered cell data"),
  make_option("--output_metada", type = 'character',
              help="Path to tsv file of filtered cell metadata")
)
opt = parse_args(OptionParser(option_list=option_list), args = commandArgs(trailingOnly = TRUE))

# Import datasets
data.counts <- read.table(
  opt$file,
  header = TRUE,
  stringsAsFactors = F,
  sep = opt$sep,
  check.names = FALSE,
  row.names = 1
)

data.counts <- data.counts[, neoplastic_ids]

QC_metrics <-
  data.frame(cell_id = colnames(data.counts),
             nGene = colSums(data.counts != 0),    # nGene : Number of detected genes for each cell
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

percentile_cutoff <- function(n, qcmetrics, variable, plot_title, ...){

# Warning I think the whole function could be refactored as described in
# https://github.com/ARTbio/mirca/pull/41#discussion_r266097289

  # Find the cutoff based on the n percentile of a distribution
  # Returns the cutoff value and plot an histogram. It requires 4 parameters :
  # n : nth percentile chosen [integer]
  # qcmetrics : matrice or dataframe that contains the number of genes and total counts metrics (QC_metrics) [object]
  # variable : either number of genes or total number of counts, from which percentile cutoff is computed [character]
  # plot_title : Title of the generated histogram (distribution of the variable) [character]
  
  var_sorted <- sort(qcmetrics[, variable]) # ascending sorting of the variable values
  percentile <- sum(var_sorted) * n / 100  # Value which covers n% of the total sum of the variable (nGene or total_counts)
  
  # Percentile loop
  new_sum <- 0
  index <- 0
  for (i in var_sorted){
    index <- index + 1
    new_sum <- new_sum + i
    if (new_sum > percentile){ # Find where your data pass the n percentile
      variable_percentile <- var_sorted[index] # the greatest value that with all other lower values
                                               # represents less than n% of the whole
      plot_hist(qcmetrics,                     # input dataframe/matrix
                variable,                      # variable used for percentile filtering
                plot_title,                    # plot title
                variable_percentile)  # Percentile threshold for filtering
      break
    }
  }
  return(variable_percentile)
}


pdf(file = opt$pdfplot)

# Determine thresholds based on percentile
if(opt$percentile_genes > 0) |(opt$percentile_counts > 0) {
  total_counts_cutoff <- percentile_cutoff(
    opt$percentile_counts,
    QC_metrics,
    "total_counts",
    "Histogram of Total counts per cell"
  )
  
  ngene_cutoff <- percentile_cutoff(
    opt$percentile_genes,
    QC_metrics,
    "nGene",
    "Histogram of Number of detected genes per cell"
  )
} else {
  ngene_cutoff <- opt$cutoff_genes
  plot_hist(QC_metrics,
            variable = "nGene",
            title = "Histogram of Number of detected genes per cell",         
            cutoff = ngene_cutoff)
  
  total_counts_cutoff <- opt$cutoff_counts
  plot_hist(QC_metrics,
            variable = "total_counts",
            title = "Histogram of Total counts per cell",
            cutoff = total_counts_cutoff)
}


# Filter out rows that didn't pass both cutoffs
QC_metrics$filtered <-
  QC_metrics$nGene <= ngene_cutoff &
  QC_metrics$total_counts <= total_counts_cutoff

# Plot the results
ggplot(QC_metrics, aes(nGene, total_counts, colour = filtered)) +
  geom_point() + scale_y_log10() +
  scale_colour_discrete(name  = "Filtered cells",
                        breaks= c(FALSE, TRUE),
                        labels= c(paste0("No (", table(QC_metrics$filtered)[1], " cells)"), 
                                  paste0("Yes (", table(QC_metrics$filtered)[2], " cells)"))) +
  xlab("Number of detected genes per cell") + ylab("Total counts per cell (log10 scale)") +
  geom_vline(xintercept = ngene_cutoff) + geom_hline(yintercept = total_counts_cutoff) +
  if(opt$percentile == T){
    ggtitle(
      paste0(
        "Filtering cells based on ",
        opt$percentile_genes,
        "th percentile of the number of genes detected by a \ncell distribution and ",
        opt$percentile_counts,
        "th percentile of the total counts per cell distribution."
      )
    )
  }else{
    ggtitle(
      paste(
        "Filtering cells that didn't detect at least",
        ngene_cutoff,
        "genes and have a minimum \nof",
        total_counts_cutoff,
        "counts"
      )
    )
  }

dev.off()

# Retrieve identifier of kept cells
kept.cells <- QC_metrics$cell_id[!QC_metrics$filtered]

data.counts <- data.counts[,kept.cells]

# Save filtered cells 
write.table(
  data.counts,
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
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
