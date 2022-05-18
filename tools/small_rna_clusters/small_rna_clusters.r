## Setup R error handling to go to stderr
options(show.error.messages = F,
        error = function() {
            cat(geterrmessage(), file = stderr()); q("no", 1, F)
        }
)
options(warn = -1)

library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(optparse)

option_list <- list(
  make_option(c("-f", "--first_dataframe"), type = "character", help = "path to first dataframe"),
  make_option("--first_plot_method", type = "character", help = "How additional data should be plotted"),
  make_option("--output_pdf", type = "character", help = "path to the pdf file with plots")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# data frames implementation

## first table
table <- read.delim(args$first_dataframe, header = T, row.names = NULL)
colnames(table)[1] <- "Dataset"
dropcol <- c("Strandness", "z.score") # not used by this Rscript and is dropped for backward compatibility
table <- table[, !(names(table) %in% dropcol)]
if (args$first_plot_method == "Counts" | args$first_plot_method == "Size") {
  table <- within(table, Counts[Polarity == "R"] <- (Counts[Polarity == "R"] * -1))
}
n_samples <- length(unique(table$Dataset))
samples <- unique(table$Dataset)
genes <- unique(table$Chromosome)
per_gene_readmap <- lapply(genes, function(x) subset(table, Chromosome == x))
per_gene_limit <- lapply(genes, function(x) c(1, unique(subset(table, Chromosome == x)$Chrom_length)))
n_genes <- length(per_gene_readmap)

## functions
plot_unit <- function(df, method = args$first_plot_method, ...) {
    p <- xyplot(Counts ~ Coordinate | factor(Dataset, levels = unique(Dataset)) + factor(Chromosome, levels = unique(Chromosome)),
               data = df,
               type = "h",
               lwd = 1.5,
               scales = list(relation = "free", x = list(rot = 0, cex = 0.7, axs = "i", tck = 0.5), y = list(tick.number = 4, rot = 90, cex = 0.7)),
               xlab = NULL, main = NULL, ylab = NULL,
               as.table = T,
               origin = 0,
               horizontal = FALSE,
               group = Polarity,
               col = c("red", "blue"),
               par.strip.text = list(cex = 0.7),
               ...)
    p <- combineLimits(p)
}

## function parameters
par_settings_firstplot <- list(layout.heights = list(top.padding = -2, bottom.padding = -2), strip.background = list(col = c("lightblue", "lightgreen")))
title_first_method <- list(Counts = "Read Counts", Coverage = "Coverage depths", Median = "Median sizes", Mean = "Mean sizes", Size = "Size Distributions")
legend_first_method <- list(Counts = "Read count", Coverage = "Coverage depth", Median = "Median size", Mean = "Mean size", Size = "Read count")
bottom_first_method <- list(Counts = "Coordinates (nucleotides)", Coverage = "Coordinates (nucleotides)", Median = "Coordinates (nucleotides)", Mean = "Coordinates (nucleotides)", Size = "Sizes of reads")

## Plotting Functions
single_plot <- function(...) {
  width <- 8.2677 * n_samples / 2
  rows_per_page <- 8
  graph_heights <- c(rep(40, 8), 10)
  pdf(file = args$output_pdf, paper = "special", height = 15, width = width)
  for (i in seq(1, n_genes, rows_per_page)) {
    start <- i
    end <- i + rows_per_page - 1
    if (end > n_genes) {
        end <- n_genes
    }
    if (end - start + 1 < 8) {
        graph_heights <- c(rep(c(40), end - start + 1), 10, rep(c(40), 8 - (end - start + 1)))
    }
    first_plot_list <- lapply(per_gene_readmap[start:end], function(x) update(useOuterStrips(plot_unit(x, par.settings = par_settings_firstplot), strip.left = strip.custom(par.strip.text = list(cex = 0.5)))))
    plot.list <- rbind(first_plot_list)
    args_list <- c(plot.list, list(nrow = rows_per_page + 1, ncol = 1, heights = unit(graph_heights, rep("mm", 9)),
                                 top = textGrob("Cluster Read Counts (Peaks in middle of clusters)", gp = gpar(cex = 1), vjust = 0, just = "top"),
                                 left = textGrob("Read Counts", gp = gpar(cex = 1), vjust = 0, hjust = 0, x = 1, y = (-0.41 / 7) * (end - start - (6.23 / 0.41)), rot = 90),
                                 sub = textGrob("Coordinates (nucleotides)", gp = gpar(cex = 1), just = "bottom", vjust = 2)
    )
    )
    do.call(grid.arrange, args_list)
  }
  devname <- dev.off()
}

# main
single_plot()
