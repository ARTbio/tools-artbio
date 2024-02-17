## Setup R error handling to go to stderr
options(
    show.error.messages = FALSE,
    error = function() {
        cat(geterrmessage(), file = stderr())
        q("no", 1, FALSE)
    }
)
warnings()
library(optparse)
library(ggplot2)
library(reshape2)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to dataframe"),
    make_option(c("-t", "--title"), type = "character", help = "Main Title"),
    make_option("--xlab", type = "character", help = "X-axis legend"),
    make_option("--ylab", type = "character", help = "Y-axis legend"),
    make_option("--sample", type = "character", help = "a space separated of sample labels"),
    make_option("--method", type = "character", help = "bedtools or pysam"),
    make_option(c("-o", "--output"), type = "character", help = "path to the pdf plot")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)
samples <- substr(args$sample, 2, nchar(args$sample) - 2)
samples <- strsplit(samples, ", ")

# data frames implementation

table <- read.delim(args$input, header = FALSE)
headers <- c("chromosome", "start", "end", "id")
for (i in seq(1, length(table) - 4)) {
    headers <- c(headers, samples[[1]][i])
    colnames(table) <- headers
}

## function
if (args$method == "bedtools") {
    cumul <- function(x, y) sum(table[, y] / (table$end - table$start) > x) / length(table$chromosome)
} else {
    cumul <- function(x, y) sum(table[, y] > x) / length(table$chromosome)
}
scale_fun <- function(x) sprintf("%.3f", x)

## end of function
## let's do a dataframe before plotting
if (args$method == "bedtools") {
    maxdepth <- trunc(max(table[, 5:length(table)] / (table$end - table$start))) + 20
} else {
    maxdepth <- trunc(max(table[, 5:length(table)])) + 20
}

graphpoints <- data.frame(1:maxdepth)
i <- 5
for (colonne in colnames(table)[5:length(colnames(table))]) {
    graphpoints <- cbind(graphpoints, mapply(cumul, 1:maxdepth, rep(i, maxdepth)))
    i <- i + 1
}
colnames(graphpoints) <- c("Depth", colnames(table)[5:length(table)])
maxfrac <- max(graphpoints[, 2:length(graphpoints)])

graphpoints <- melt(graphpoints, id.vars = "Depth", variable.name = "Samples", value.name = "sample_value")

## GRAPHS

pdf(file = args$output)
ggplot(data = graphpoints, aes(x = Depth, y = sample_value, colour = Samples)) +
    geom_line(size = 1) +
    scale_x_continuous(trans = "log2", breaks = 2^(seq(0, log(maxdepth, 2)))) +
    scale_y_continuous(breaks = seq(0, maxfrac, by = maxfrac / 10), labels = scale_fun) +
    labs(x = args$xlab, y = args$ylab, title = args$title) +
    theme(
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(colour = "blue", size = 7)
    )

devname <- dev.off()
