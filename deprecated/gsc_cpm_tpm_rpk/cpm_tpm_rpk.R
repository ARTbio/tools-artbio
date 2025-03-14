if (length(commandArgs(TRUE)) == 0) {
  system("Rscript cpm_tpm_rpk.R -h", intern = F)
  q("no")
}


# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)

#Arguments
option_list = list(
  make_option(
    c("-d", "--data"),
    default = NA,
    type = 'character',
    help = "Input file that contains count values to transform"
  ),
  make_option(
    c("-t", "--type"),
    default = 'cpm',
    type = 'character',
    help = "Transformation type, either cpm, tpm, rpk or none[default : '%default' ]"
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
    help = "Consider first line as header ? [default : '%default' ]"
  ),
  make_option(
    c("-f", "--gene"),
    default = NA,
    type = 'character',
    help = "Two column of gene length file"
  ),
  make_option(
    c("-a", "--gene_sep"),
    default = '\t',
    type = 'character',
    help = "Gene length file separator [default : '%default' ]"
  ),
  make_option(
    c("-b", "--gene_header"),
    default = TRUE,
    type = 'logical',
    help = "Consider first line of gene length as header ? [default : '%default' ]"
  ),
  make_option(
    c("-l", "--log"),
    default = FALSE,
    type = 'logical',
    help = "Should be log transformed as well ? (log2(data +1)) [default : '%default' ]"
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
if (opt$gene_sep == "tab") {opt$gene_sep = "\t"}

cpm <- function(count) {
  t(t(count) / colSums(count)) * 1000000
}


rpk <- function(count, length) {
  count / (length / 1000)
}

tpm <- function(count, length) {
  RPK = rpk(count, length)
  perMillion_factor = colSums(RPK) / 1000000
  TPM = RPK / perMillion_factor
  return(TPM)
}

data = read.table(
  opt$data,
  check.names = FALSE,
  header = opt$colnames,
  row.names = 1,
  sep = opt$sep
)

if (opt$type == "tpm" | opt$type == "rpk") {
  gene_length = as.data.frame(
    read.table(
      opt$gene,
      h = opt$gene_header,
      row.names = 1,
      sep = opt$gene_sep
    )
  )
  gene_length = as.data.frame(gene_length[match(rownames(data), rownames(gene_length)), ], rownames(data))
}


if (opt$type == "cpm")
  res = cpm(data)
if (opt$type == "tpm")
  res = as.data.frame(apply(data, 2, tpm, length = gene_length), row.names = rownames(data))
if (opt$type == "rpk")
  res = as.data.frame(apply(data, 2, rpk, length = gene_length), row.names = rownames(data))
if (opt$type == "none")
  res = data
colnames(res) = colnames(data)


if (opt$log == TRUE) {
  res = log2(res + 1)
}

write.table(
  cbind(Genes = rownames(res), res),
  opt$out,
  col.names = opt$colnames,
  row.names = F,
  quote = F,
  sep = "\t"
)












