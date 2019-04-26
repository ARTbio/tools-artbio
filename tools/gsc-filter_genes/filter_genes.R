#########################
#      filter genes     #
#########################

# Third step of the signature-based workflow
# Filter low expressed genes

#Example of command (used for generate output file) :
#Rscript 3-filter_genes.R -d ../2-log2CPM1P/log2CPM1p.tsv -o .

# Load necessary packages (install them if it's not the case)
requiredPackages = c('optparse')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE, quietly = T)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(suppressMessages(library(p, character.only = TRUE)))
}


#Arguments
option_list = list(
  make_option(
    c("-d", "--data"),
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
    help = "Consider first line as header ? [default : '%default' ]"
  ),
  make_option(
    c("-p", "--percent"),
    default = 0.03,
    type = 'numeric',
    help = "Include genes with detected expression in at least \
    this fraction of cells [default : '%default' ]"
  ),
  # make_option(
  #   c("-m", "--min.cells"),
  #   default = 3,
  #   type = 'integer',
  #   help = "Include genes with detected expression in at least n cells [default : '%default' ]"
  # ),
  make_option(
    c("-o", "--out"),
    default = "~",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)


opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$data == "" & !(opt$help)) {
  stop("At least one argument must be supplied (count data --data option).\n",
       call. = FALSE)
}

#Open files
data.counts <- read.table(
  opt$data,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

#Search for genes that are expressed in a certain percent of cells
kept_genes <- rowSums(data.counts != 0) > (opt$percent * ncol(data.counts))
# #Search for genes that are expressed in a certain number of cells
# kept_genes <- rowSums(data.counts != 0) >= opt$min.cells

#Filter matrix
data.counts <- data.counts[kept_genes,]

#Save filtered matrix
write.table(
  data.counts,
  paste(opt$out, "filterGenes.tsv", sep = "/"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)