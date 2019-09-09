# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()

library(optparse)
library(scran)
print("ok")

#Arguments
option_list = list(
  make_option(
    c("-d", "--data"),
    default = NA,
    type = 'character',
    help = "Input file that contains count values to transform"
  ),
  make_option(
    c("-s", "--sep"),
    default = '\t',
    type = 'character',
    help = "File separator [default : '%default' ]"
  ),
  make_option(
    c("-o", "--out"),
    default = "res.tab",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
print(opt$sep)
data = read.table(
  opt$data,
  check.names = FALSE,
  header = TRUE,
  row.names = 1,
  sep = opt$sep
)

## Import data as a SingleCellExperiment object
sce <- SingleCellExperiment(list(counts=as.matrix(data)))

## Compute sum factors
sce <- computeSumFactors(sce)

sce <- normalize(sce)

logcounts <- data.frame(genes = rownames(sce), logcounts(sce), check.names = F)


write.table(
  logcounts,
  opt$out,
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)
