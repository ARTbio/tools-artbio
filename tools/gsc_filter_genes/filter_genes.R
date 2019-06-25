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
    "--percentile_detection",
    default = 0,
    type = 'numeric',
    help = "Include genes with detected expression in at least \
    this fraction of cells [default : '%default' ]"
  ),
  make_option(
    "--absolute_detection",
    default = 0,
    type = 'numeric',
    help = "Include genes with detected expression in at least \
    this number of cells [default : '%default' ]"
  ),
  make_option(
    c("-o", "--output"),
    default = NA,
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))
if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

# Open files
data.counts <- read.table(
  opt$input,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

# note the [if else] below, to handle percentile_detection=absolute_detection=0
# Search for genes that are expressed in a certain percent of cells
if (opt$percentile_detection > 0) {
kept_genes <- rowSums(data.counts != 0) >= (opt$percentile_detection * ncol(data.counts))
} else {

# Search for genes that are expressed in more than an absolute number of cells
kept_genes <- rowSums(data.counts != 0) >= (opt$absolute_detection)
}

# Filter matrix
data.counts <- data.counts[kept_genes,]
data.counts <- cbind(Genes=rownames(data.counts), data.counts)

# Save filtered matrix
write.table(
  data.counts,
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)