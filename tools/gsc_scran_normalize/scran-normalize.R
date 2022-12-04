# load packages that are provided in the conda env
options(show.error.messages=F,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
    }
)
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()

library(optparse)
library(scran)

# Arguments
option_list <- list(
  make_option(
    c("-d", "--data"),
    default = NA,
    type = "character",
    help = "Input file that contains count values to transform"
  ),
  make_option(
    c("-s", "--sep"),
    default = "\t",
    type = "character",
    help = "File separator [default : '%default' ]"
  ),
  make_option(
    "--cluster",
    default = FALSE,
    action = "store_true",
    type = "logical",
    help = "Whether to calculate the size factor per cluster or on all cell"
  ),
  make_option(
    c("-m", "--method"),
    default = "hclust",
    type = "character",
    help = "The clustering method to use for grouping cells into cluster : hclust or igraph [default : '%default' ]"
  ),
  make_option(
    "--size",
    default = 100,
    type = "integer",
    help = "Minimal number of cells in each cluster : hclust or igraph [default : '%default' ]"
  ),
  make_option(
    c("-o", "--out"),
    default = "res.tab",
    type = "character",
    help = "Output name [default : '%default' ]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {
  opt$sep = "\t"
  }

data <- read.table(
  opt$data,
  check.names = FALSE,
  header = TRUE,
  row.names = 1,
  sep = opt$sep
)

## Import data as a SingleCellExperiment object
sce <- SingleCellExperiment(list(counts = as.matrix(data)))

if (opt$cluster) {
  clusters <- quickCluster(sce, min.size = opt$size, method = opt$method)

  ## Compute sum factors
  sce <- computeSumFactors(sce, cluster = clusters)
} else {

  ## Compute sum factors
  sce <- computeSumFactors(sce)
}

sce <- normalize(sce)

logcounts <- data.frame(genes = rownames(sce), round(logcounts(sce), digits = 5), check.names = FALSE)


write.table(
  logcounts,
  opt$out,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
