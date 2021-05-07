## Setup R error handling to go to stderr
options( show.error.messages=F,
         error = function () {
           cat( geterrmessage(), file = stderr() ); q( "no", 1, F )
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
  make_option(c("-i", "--ymin"), type = "double", help = "set min ylimit. e.g. '-100.0'"),
  make_option(c("-a", "--ymax"), type = "double", help = "set max ylimit. e.g. '100.0'"),
  make_option(c("-f", "--first_dataframe"), type = "character", help = "path to first dataframe"),
  make_option(c("-e", "--extra_dataframe"), type = "character", help="path to additional dataframe"),
  make_option(c("-n", "--normalization"), type = "character", help = "space-separated normalization/size factors"),
  make_option("--first_plot_method", type = "character", help = "How additional data should be plotted"),
  make_option("--extra_plot_method", type = "character", help = "How additional data should be plotted"),
  make_option("--global", type = "character", help = "data should be plotted as global size distribution"),
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
  table <- within(table, Counts[Polarity=="R"] <- (Counts[Polarity=="R"]*-1))
}
n_samples <- length(unique(table$Dataset))
samples <- unique(table$Dataset)
if (args$normalization != "") {
  norm_factors <- as.numeric(unlist(strsplit(args$normalization, " ")))
} else {
  norm_factors <- rep(1, n_samples)
}
if (args$first_plot_method == "Counts" | args$first_plot_method == "Size" | args$first_plot_method == "Coverage") {
  i <- 1
  for (sample in samples) {
    # Warning Here the column is hard coded as the last column (dangerous)
    # because its name changes with the method
    table[, length(table)][table$Dataset == sample] <- table[, length(table)][table$Dataset == sample] * norm_factors[i]
    i <- i + 1
  }
}
genes <- unique(table$Chromosome)
per_gene_readmap <- lapply(genes, function(x) subset(table, Chromosome == x))
per_gene_limit <- lapply(genes, function(x) c(1, unique(subset(table, Chromosome == x)$Chrom_length)) )
n_genes <- length(per_gene_readmap)

# second table
if (args$extra_plot_method != "") {
  extra_table <- read.delim(args$extra_dataframe, header = T, row.names = NULL)
  colnames(extra_table)[1] <- "Dataset"
  dropcol <- c("Strandness", "z.score")
  table <- table[, !(names(table) %in% dropcol)]
  if (args$extra_plot_method == "Counts" | args$extra_plot_method == "Size") {
    extra_table <- within(extra_table, Counts[Polarity == "R"] <- (Counts[Polarity == "R"] * -1))
  }
  if (args$extra_plot_method == "Counts" | args$extra_plot_method == "Size" | args$extra_plot_method == "Coverage") {
    i <- 1
    for (sample in samples) {
      extra_table[, length(extra_table)][extra_table$Dataset == sample] <- extra_table[, length(extra_table)][extra_table$Dataset == sample] * norm_factors[i]
      i = i + 1
    }
  }
  per_gene_size = lapply(genes, function(x) subset(extra_table, Chromosome == x))
}

## functions
globalbc = function(df, global = "", ...) {
  if (global == "yes") {
    bc <- barchart(Counts ~ as.factor(Size) | factor(Dataset, levels = unique(Dataset)),
                   data = df, origin = 0,
                   horizontal = FALSE,
                   col = c("darkblue"),
                   scales = list(y = list(tick.number = 4, rot = 90, relation = "same", cex = 0.5, alternating = T), x = list(rot = 0, cex = 0.6, tck = 0.5, alternating = c(3, 3))),
                   xlab = list(label = bottom_first_method[[args$first_plot_method]], cex = .85),
                   ylab = list(label=legend_first_method[[args$first_plot_method]], cex = .85),
                   main = title_first_method[[args$first_plot_method]],
                   layout = c(2, 6), newpage = T,
                   as.table = TRUE,
                   aspect = 0.5,
                   strip = strip.custom(par.strip.text = list(cex = 1), which.given = 1, bg = "lightblue"),
                   ...
    )
  } else {
    bc <- barchart(Counts ~ as.factor(Size) | factor(Dataset, levels = unique(Dataset)),
                   data = df, origin = 0,
                   horizontal = FALSE,
                   group = Polarity,
                   stack = TRUE,
                   col = c('red', 'blue'),
                   scales = list(y = list(tick.number = 4, rot= 90, relation = "same", cex = 0.5, alternating = T), x = list(rot = 0, cex = 0.6, tck = 0.5, alternating = c(3,3))),
                   xlab = list(label = bottom_first_method[[args$first_plot_method]], cex = .85),
                   ylab = list(label = legend_first_method[[args$first_plot_method]], cex = .85),
                   main = title_first_method[[args$first_plot_method]],
                   layout = c(2, 6), newpage = T,
                   as.table = TRUE,
                   aspect = 0.5,
                   strip = strip.custom(par.strip.text = list(cex = 1), which.given = 1, bg = "lightblue"),
                   ...
    )
  }
  return(bc)
}
plot_unit <- function(df, method = args$first_plot_method, ...) {
  if (exists("ymin", where = args)) {
    min <- args$ymin
  } else {
    min <- ""
  }
  if ((exists("ymax", where = args))) {
    max <- args$ymax
  } else {
    max <- ""
  }
  ylimits = c(min, max)
  if (method == 'Counts') {
    p <- xyplot(Counts ~ Coordinate | factor(Dataset, levels = unique(Dataset)) + factor(Chromosome, levels = unique(Chromosome)),
               data = df,
               type = 'h',
               lwd = 1.5,
               scales = list(relation="free", x = list(rot = 0, cex = 0.7, axs = "i", tck = 0.5), y = list(tick.number = 4, rot = 90, cex = 0.7)),
               xlab = NULL, main = NULL, ylab = NULL, ylim = ylimits,
               as.table = T,
               origin = 0,
               horizontal = FALSE,
               group = Polarity,
               col = c("red", "blue"),
               par.strip.text = list(cex = 0.7),
               ...)
    p <- combineLimits(p)
  } else if (method != "Size") {
    p <- xyplot(eval(as.name(method)) ~ Coordinate | factor(Dataset, levels = unique(Dataset)) + factor(Chromosome, levels = unique(Chromosome)),
               data = df,
               type = ifelse(method =='Coverage', 'l', 'p'),
               pch = 19,
               cex = 0.35,
               scales = list(relation="free", x = list(rot = 0, cex = 0.7, axs = "i", tck = 0.5), y = list(tick.number = 4, rot = 90, cex = 0.7)),
               xlab = NULL, main = NULL, ylab = NULL, ylim = ylimits,
               as.table = T,
               origin = 0,
               horizontal = FALSE,
               group = Polarity,
               col = c("red", "blue"),
               par.strip.text = list(cex = 0.7),
               ...)
    p <- combineLimits(p)
  } else {
    p <- barchart(Counts ~ as.factor(Size) | factor(Dataset, levels = unique(Dataset)) + Chromosome, data = df, origin = 0,
                 horizontal = FALSE,
                 group = Polarity,
                 stack = TRUE,
                 col = c('red', 'blue'),
                 scales = list(y = list(rot = 90, relation = "free", cex = 0.7), x = list(rot = 0, cex = 0.7, axs = "i", tck=c(1, 0))),
                 xlab = NULL,
                 ylab = NULL,
                 main = NULL,
                 as.table = TRUE,
                 par.strip.text = list(cex = 0.6),
                 ...)
    p <- combineLimits(p)
  }
  return(p)
}


## function parameters

par.settings.firstplot = list(layout.heights = list(top.padding = -2, bottom.padding = -2), strip.background = list(col = c("lightblue", "lightgreen")))
par.settings.secondplot = list(layout.heights = list(top.padding = -1, bottom.padding = -1), strip.background = list(col = c("lightblue", "lightgreen")))
title_first_method = list(Counts = "Read Counts", Coverage = "Coverage depths", Median = "Median sizes", Mean = "Mean sizes", Size = "Size Distributions")
title_extra_method = list(Counts = "Read Counts", Coverage = "Coverage depths", Median = "Median sizes", Mean = "Mean sizes", Size = "Size Distributions")
legend_first_method = list(Counts= "Read count", Coverage = "Coverage depth", Median = "Median size", Mean = "Mean size", Size = "Read count")
legend_extra_method = list(Counts = "Read count", Coverage = "Coverage depth", Median = "Median size", Mean = "Mean size", Size = "Read count")
bottom_first_method = list(Counts = "Coordinates (nucleotides)", Coverage = "Coordinates (nucleotides)", Median = "Coordinates (nucleotides)", Mean = "Coordinates (nucleotides)", Size = "Sizes of reads")
bottom_extra_method = list(Counts = "Coordinates (nucleotides)", Coverage = "Coordinates (nucleotides)", Median = "Coordinates (nucleotides)", Mean = "Coordinates (nucleotides)", Size = "Sizes of reads")

## Plotting Functions

double_plot <- function(...) {
  page_height <- 15
  rows_per_page <- 10
  graph_heights <- c(40, 30, 40, 30, 40, 30, 40, 30, 40, 30, 10)
  page_width <- 8.2677 * n_samples / 2
  pdf(file = args$output_pdf, paper = "special", height = page_height, width = page_width)
  for (i in seq(1, n_genes, rows_per_page / 2)) {
    start <- i
    end <- i + rows_per_page / 2 - 1
    if (end > n_genes) {
      end <- n_genes
      }
    if (end-start+1 < 5) {
      graph_heights <- c(rep(c(40,30), end - start + 1), 10, rep(c(40, 30), 5 - (end - start + 1)))}
    first_plot_list <- lapply(per_gene_readmap[start:end], function(x) update(useOuterStrips(plot_unit(x, par.settings = par.settings.secondplot), strip.left = strip.custom(par.strip.text = list(cex = 0.5)))))
    second_plot_list <- lapply(per_gene_size[start:end], function(x) update(useOuterStrips(plot_unit(x, method = args$extra_plot_method, par.settings = par.settings.firstplot), strip.left = strip.custom(par.strip.text = list(cex = 0.5)), strip = FALSE)))
    plot.list <- rbind(first_plot_list, second_plot_list)
    args_list <- c(plot.list, list(nrow = rows_per_page + 1, ncol = 1, heights = unit(graph_heights, rep("mm", 11)),
                                 top = textGrob(paste(title_first_method[[args$first_plot_method]], "and", title_extra_method[[args$extra_plot_method]]), gp = gpar(cex = 1), vjust = 0, just = "top"),
                                 left = textGrob(paste(legend_first_method[[args$first_plot_method]], "/", legend_extra_method[[args$extra_plot_method]]), gp = gpar(cex = 1), vjust = 0, hjust = 0, x = 1, y = (-0.38 / 4) * (end - start - (3.28 / 0.38)), rot = 90),
                                 sub = textGrob(paste(bottom_first_method[[args$first_plot_method]], "/", bottom_extra_method[[args$extra_plot_method]]), gp = gpar(cex=1), just = "bottom", vjust = 2)
    )
    )
    do.call(grid.arrange, args_list)
  }
  devname <- dev.off()
}


single_plot <- function(...) {
  width <- 8.2677 * n_samples / 2
  rows_per_page <- 8
  graph_heights <- c(rep(40, 8), 10)
  pdf(file = args$output_pdf, paper = "special", height = 15, width = width)
  for (i in seq(1, n_genes, rows_per_page)) {
    start = i
    end <- i + rows_per_page - 1
    if (end > n_genes) {
      end = n_genes
      }
    if (end - start + 1 < 8) {
      graph_heights <- c(rep(c(40), end - start + 1), 10, rep(c(40), 8 - (end - start + 1)))}
    first_plot_list <- lapply(per_gene_readmap[start:end], function(x) update(useOuterStrips(plot_unit(x, par.settings = par.settings.firstplot), strip.left = strip.custom(par.strip.text = list(cex = 0.5)))))
    plot.list <- rbind(first_plot_list)
    args_list <- c(plot.list, list(nrow = rows_per_page + 1, ncol = 1, heights = unit(graph_heights, rep("mm", 9)),
                                 top = textGrob(title_first_method[[args$first_plot_method]], gp=gpar(cex = 1), vjust = 0, just = "top"),
                                 left = textGrob(legend_first_method[[args$first_plot_method]], gp = gpar(cex = 1), vjust = 0, hjust = 0, x = 1, y = (-0.41/7) * (end - start - (6.23/0.41)), rot = 90),
                                 sub = textGrob(bottom_first_method[[args$first_plot_method]], gp = gpar(cex = 1), just = "bottom", vjust = 2)
    )
    )
    do.call(grid.arrange, args_list)
  }
  devname <- dev.off()
}

# main

if (args$extra_plot_method != "") {
  double_plot() }
if (args$extra_plot_method == "" & !exists("global", where = args)) {
  single_plot()
}
if (exists("global", where = args)) {
  pdf(file = args$output, paper = "special", height = 11.69)
  table <- within(table, Counts[Polarity == "R"] <- abs(Counts[Polarity = ="R"]))
  library(reshape2)
  ml <- melt(table, id.vars = c("Dataset", "Chromosome", "Polarity", "Size"))
  if (args$global == "nomerge") {
    castml <- dcast(ml, Dataset + Polarity+Size ~ variable, function(x) sum(x))
    castml <- within(castml, Counts[Polarity == "R"] <- (Counts[Polarity == "R"] * -1))
    bc = globalbc(castml, global="no")
  } else {
    castml <- dcast(ml, Dataset+Size ~ variable, function(x) sum(x))
    bc <- globalbc(castml, global = "yes")
  }
  plot(bc)
  devname <- dev.off()
}
