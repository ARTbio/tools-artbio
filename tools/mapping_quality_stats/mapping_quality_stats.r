## Setup R error handling to go to stderr
options(show.error.messages = F,
            error = function() {
                cat(geterrmessage(), file = stderr()); q("no", 1, F)
            }
       )
warnings()
library(optparse)
library(ggplot2)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to tabular file"),
    make_option(c("-o", "--output"), type = "character", help = "path to the pdf plot")
    )
 
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# data frame implementation
table <- read.delim(args$input, header = T)
colnames(table) <- c("MAPQ", "Counts")


# Barplot
pdf(file = args$output)
ggplot(data, aes(x=Counts, y=MAPQ)) + 
    geom_bar(stat = "identity")
devname <- dev.off()
