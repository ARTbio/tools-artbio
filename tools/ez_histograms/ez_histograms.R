library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(vtable)
library(optparse)

options(show.error.messages = FALSE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)

loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()

option_list <- list(
  make_option(
    c("-f", "--file"),
    default = NA,
    type = "character",
    help = "Input file that contains count values to transform"
  ),
  make_option(
    c("-d", "--profile"),
    default = "count",
    type = "character",
    help = "Whether y-axis shows absolute counts or density: 'count' or 'density' [default : '%default' ]"
  ),
  make_option(
    "--xscale",
    default = "cartesian",
    type = "character",
    help = "Whether x-axis is 'cartesian', 'log2' or 'log10' [default : '%default' ]"
  ),
  make_option(
    "--yscale",
    default = "cartesian",
    type = "character",
    help = "Whether y-axis is 'cartesian', 'log2' or 'log10' [default : '%default' ]"
  ),
  make_option(
    c("-p", "--pdf"),
    default = "histograms.pdf",
    type = "character",
    help = "Output pdf file name [default : '%default' ]"
  ),
  make_option(
    c("-s", "--summary"),
    default = "summary.tsv",
    type = "character",
    help = "statistics summary file name [default : '%default' ]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                  args = commandArgs(trailingOnly = TRUE))

plot_histograms <- function(mdata, profile = "count", xscale = "cartesian", yscale = "cartesian", bins = 30) {
  if (profile == "count") {
    # count histogram
    p <- ggplot(mdata, aes(x = value, fill = variable, color = variable, y = after_stat(count)), show.legend = FALSE) +
      geom_histogram(bins = bins) + theme(legend.position = "none")
    if (xscale == "cartesian") {
      if (yscale == "log2") {
        p <- p + scale_y_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
      } else {
        if (yscale == "log10") {
          p <- p + scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))
        }
      }
    }
    if (xscale == "log2") {
      p <- p + scale_x_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
      if (yscale == "log2") {
        p <- p + scale_y_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
      } else {
        if (yscale == "log10") {
          p <- p + scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))
        }
      }
    }
    if (xscale == "log10") {
      p <- p + scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))
      if (yscale == "log2") {
        p <- p + scale_y_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
      } else {
        if (yscale == "log10") {
          p <- p + scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))
        }
      }
    }
  }

  if (profile == "density") {
    # density histogram
    p <- ggplot(mdata, aes(x = value, fill = variable, color = variable)) +
      geom_density() + theme(legend.position = "none")
    if (xscale == "log2") {
      p <- p + scale_x_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
    }
    if (xscale == "log10") {
      p <- p + scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))
    }
  }
  return(p)
}

test_header <- function(file) {
  data <- read.delim(file = file, header = FALSE, row.names = 1, nrows = 2)
  if(all(is.na(as.numeric(data[1, seq_len(ncol(data))])))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

test_rownames <- function(file) {
  data <- read.delim(file = file, header = FALSE, row.names = NULL, nrows = 2)
  if (is.na(as.numeric(data[2, 1]))) {
    return(1)
  } else {
    return(NULL)
  }
}

##### prepare input data
data <- read.delim(file = opt$file, header = test_header(opt$file), row.names = test_rownames(opt$file))
data <- data %>% select(where(is.numeric))  # remove non numeric columns
mdata <- melt(data)

##### main

# determine optimal number of bins (Sturges’ Rule)
bins <- ceiling(log2(nrow(data)) + 1)
# plot
p <- plot_histograms(mdata, profile = opt$profile, xscale = opt$xscale, bins = bins, yscale = opt$yscale)

# determine optimal width for the graph
width <- length(data)
width <- case_when(
  width == 1 ~ 14 / 3,
  width == 2 ~ (2 / 3) * 14,
  TRUE ~ 14
)
# determine optimal height for the graph
height <- length(data)
height <- case_when(
  height <= 3 ~ 3,
  height <= 6 ~ 6,
  TRUE ~ (floor(height / 3) + 1) * 3
)
# determine optimal number of col for the graph
ncol <- length(data)
ncol <- case_when(
  ncol == 1 ~ 1,
  ncol == 2 ~ 2,
  TRUE ~ 3
)
pdf(opt$pdf, width = width, height = height)
print(p + facet_wrap(~variable, ncol = ncol, scales = "free"))
dev.off()

# Summary statistics with vtable package
summary_df <- sumtable(data, digits = 8, out = "return", add.median = TRUE,
                       summ.names = c("N", "Mean", "Std. Dev.", "Min", "Pctl. 25",
                                      "Median", "Pctl. 75", "Max"))
write.table(summary_df, file = opt$summary, sep = "\t", quote = FALSE, row.names = FALSE)
