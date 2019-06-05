# Performs multi-correlation analysis between the vectors of gene expressions
# in single cell RNAseq libraries and the vectors of signature scores in these
# same single cell RNAseq libraries.

# Example of command
# Rscript correlations_with_signature.R --expression_file <expression_data.tsv>
#                                       -m <cell_signature_scores.tsv>
#                                       -g ../5-differential_analysis/diffAnalysisGeneMetadata.tsv
#                                       --score Signature_score -o .

# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
requiredPackages = c('optparse', 'Hmisc', "heatmaply")
warnings()
library(optparse)
# both next packages are in conda-forge (r-hmisc and r-heatmaply)
library(Hmisc)
library(heatmaply)

#Arguments
option_list = list(
  make_option(
    "--sep",
    default = '\t',
    type = 'character',
    help = "File separator, must be the same for all input files [default : '%default' ]"
  ),
  make_option(
    c("-c", "--colnames"),
    default = TRUE,
    type = 'logical',
    help = "Consider first lines as header (must stand for all input files)? [default : '%default' ]"
  ),  
  make_option(
    --expression_file,
    default = NA,
    type = 'character',
    help = "Input file that contains log2(CPM +1) expression values"
  ),
  make_option(
    --signatures_file,
    default = NA,
    type = 'character',
    help = "Input file that contains cell signature"
  ),
  make_option(
    "--score",
    default = "score",
    type = 'character',
    help = "Column name of signature score in cell signatures file [default : '%default' ]"
  ), 
  make_option(
    "--DE_genes",
    default = NA,
    type = 'character',
    help = "Input file that contains genes metadata"
  ),
  make_option(
    c("-a", "--autocorr"),
    default = FALSE,
    type = 'logical',
    help = "Output the gene-gene expression correlation matrices [default : '%default' ]"
  ),
  make_option(
    --output,
    default = "output",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

# if (opt$file == "" | opt$metadata == "" & !(opt$help)) {
#   stop("At least tw arguments must be supplied (count data --file option and cell metadata --metadata option).\n",
#        call. = FALSE)
# }

# Open files
data.counts <- read.table(
  opt$expression_file,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)
metadata <- read.delim(
  opt$signatures_file,
  header = T,
  stringsAsFactors = F,
  sep = opt$sep,
  check.names = F,
  row.names = 1
)
gene_metadata <- read.delim(
  opt$DE_genes,
  header = T,
  stringsAsFactors = F,
  sep = opt$sep,
  check.names = F,
  row.names = 1
)

# Retrieve significantly differentially expressed genes
DE_genes <- rownames(subset(gene_metadata, Significant == TRUE)) # here a header name is imposed, to do :  select on padj !

# Filter expression matrix
data <- data.counts[DE_genes,]

# Add signature score to expression matrix
data <- rbind(t(subset(metadata, select = c(opt$score))), data[,rownames(metadata)])

# Gene correlation
gene_corr <- rcorr(t(data), type = "pearson") # transpose because we correlate genes, not cells

#Gene correlation with signature score
gene_corr_score <- cbind.data.frame(r = gene_corr$r[, opt$score], 
                                    P = gene_corr$P[, opt$score])

#Heatmap
heatmaply(
  subset(gene_corr_score, select = r),
  xlab = "Signature Score", 
  ylab = "Differentially Expressed Genes",
  main = paste0("Correlation between Signature and significantly DE genes"),
  dendrogram = "row",
  colors = RdBu,
  showticklabels = c(F,F),
  limits = c(-1,1),
  margins = c(40,NA,NA,NA),
  file = paste0(opt$out, "/CorrelationHeatmap.html")
)


#Save files
write.table(
  gene_corr_score,
  paste0(opt$out, "/GeneCorrelationWithSignature.tsv"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)

if(opt$autocorr == T){
  write.table(
    gene_corr$r[-1,-1],
    paste0(opt$out, "/GeneGeneCorrelation.tsv"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = T
  )
  write.table(
    gene_corr$P[-1,-1],
    paste0(opt$out, "/GeneGeneCorrelationPvalues.tsv"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = T
  )
}
