# ########################
#      filter genes     #
# ########################

# Filter out low expressed genes

# Example of command (used for generate output file) :
# Rscript filter_genes.R -f <input file> -o <output file>

# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
library(optparse)

# Arguments
option_list = list(
  make_option(
    c("-f", "--input"),
    default = NA,
    type = 'character',
    help = "Input file that contains count values to filter"
  ),
  make_option(
    c("-s", "--sep"),
    default = '\t',
    type = 'character',
    help = "File separator [default : '%default' ]"
  ),
  make_option(
    c("-c", "--colnames"),
    default = TRUE,
    type = 'logical',
    help = "first line is a header [default : '%default' ]"
  ),
  make_option(
    "--detection",
    default = 0.05,
    type = 'numeric',
    help = "Include genes with detected expression in at least \
    this fraction of cells [default : '%default' ]"
  ),
  make_option(
    c("-o", "--output"),
    default = NONE
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

# Open files
data.counts <- read.table(
  opt$input,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

# Search for genes that are expressed in a certain percent of cells
kept_genes <- rowSums(data.counts != 0) > (opt$detection * ncol(data.counts))

# Filter matrix
data.counts <- data.counts[kept_genes,]

# Save filtered matrix
write.table(
  data.counts,
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)