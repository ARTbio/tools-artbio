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
library(ggplot2)
library(reshape2)
library(Rtsne)



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
    help = "Transformation type, either cpm, tpm or rpk [default : '%default' ]"
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
  ),
  make_option(
    "--tsne",
    default = FALSE,
    type = 'logical',
    help = "performs T-SNE [default : '%default' ]"
  ),
  make_option(
    "--seed",
    default = 42,
    type = 'integer',
    help = "Seed value for reproducibility [default : '%default' ]"
  ),
  make_option(
    "--perp",
    default = 5.0,
    type = 'numeric',
    help = "perplexity [default : '%default' ]"
  ),
  make_option(
    "--theta",
    default = 1.0,
    type = 'numeric',
    help = "theta [default : '%default' ]"
  ),
  make_option(
    "--legend",
    default = TRUE,
    type = 'logical',
    help = "Legend options [default : '%default' ]"
  ),
  make_option(
    c("-D", "--tsne_out"),
    default = "tsne.pdf",
    type = 'character',
    help = "T-SNE pdf [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$data == "" & !(opt$help)) {
  stop("At least one argument must be supplied (count data).\n",
       call. = FALSE)
} else if ((opt$type == "tpm" | opt$type == "rpk") & opt$gene == "") {
  stop("At least two arguments must be supplied (count data and gene length file).\n",
       call. = FALSE)
} else if (opt$type != "tpm" & opt$type != "rpk" & opt$type != "cpm") {
  stop("Wrong transformation requested (--type option) must be : cpm, tpm or rpk.\n",
       call. = FALSE)
}

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
colnames(res) = colnames(data)


if (opt$log == TRUE) {
  res = log2(res + 1)
}

write.table(
  cbind(Features = rownames(res), res),
  opt$out,
  col.names = opt$colnames,
  row.names = F,
  quote = F,
  sep = "\t"
)

## 
if (opt$tsne == TRUE) {
  df = cpm(data)
  # filter and transpose df for tsne
#  rownames(df) = df[,1]
#  df = df[,-1] # we remove the first item column for the transposition
  df = df[rowSums(df) != 0,] # remove lines without information (with only 0 counts)
#  first_column_name = colnames(df)[1] # save the column name for latter
  tdf = t(df)
  tdf =  log2(tdf + 1)
  # make tsne and plot results
## Show/hide legend
  if (opt$legend == TRUE) {
    gg_legend = NULL } else {
    gg_legend = theme(legend.position="none")
  }
  gg_legend = theme(legend.position="none")
  set.seed(opt$seed) ## Sets seed for reproducibility
  # Run TSNE
  tsne_out <- Rtsne(tdf, perplexity=opt$perp, theta=opt$theta) # 
  embedding <- as.data.frame(tsne_out$Y)
  embedding$Class <- as.factor(sub("Class_", "", rownames(tdf)))
  ggplot(embedding, aes(x=V1, y=V2, color=Class)) +
    geom_point(size=1.25) +
    geom_text(aes(label=Class),hjust=-0.2, vjust=-0.5, size=2.5) +
    gg_legend +
    xlab("") +
    ylab("") +
    ggtitle('t-SNE of data (log2CPM transformed)')
  ggsave(file=opt$tsne_out, device="pdf")
}
  











