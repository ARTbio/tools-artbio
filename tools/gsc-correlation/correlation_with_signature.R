####################
#    Correlation   #
#     analysis     #
####################

# Sixth step of the signature-based workflow
# Perform a correlation analysis between the
# signature score and DEG expression.

# Example of command
# Rscript 6-correlation.R -f ../3-filter_genes/filterGenes.tsv -m ../4-signature_score/signatureScoreMetadata.tsv -g ../5-differential_analysis/diffAnalysisGeneMetadata.tsv --score Signature_score -o .

# Load necessary packages (install them if it's not the case)
requiredPackages = c('optparse', 'Hmisc', "heatmaply")
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE, quietly = T)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(suppressMessages(library(p, character.only = TRUE)))
}

#Arguments
option_list = list(
  make_option(
    c("-f", "--file"),
    default = NA,
    type = 'character',
    help = "Input file that contains log2(CPM +1) values"
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
    c("-m", "--metadata"),
    default = NA,
    type = 'character',
    help = "Input file that contains cells metadata"
  ),
  make_option(
    c("-d", "--delimiter"),
    default = "\t",
    type = 'character',
    help = "Column separator for metadata file [default : '%default' ]"
  ),
  make_option(
    "--score",
    default = "score",
    type = 'character',
    help = "Column name of signature score in gene metadata file [default : '%default' ]"
  ), 
  make_option(
    c("-g", "--gene_metadata"),
    default = NA,
    type = 'character',
    help = "Input file that contains genes metadata"
  ),
  make_option(
    "--gene_delimiter",
    default = "\t",
    type = 'character',
    help = "Column separator for gene metadata file [default : '%default' ]"
  ),
  make_option(
    c("-a", "--autocorr"),
    default = FALSE,
    type = 'logical',
    help = "Output the gene-gene expression correlation matrices [default : '%default' ]"
  ),
  make_option(
    c("-o", "--out"),
    default = "~",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$file == "" | opt$metadata == "" & !(opt$help)) {
  stop("At least tw arguments must be supplied (count data --file option and cell metadata --metadata option).\n",
       call. = FALSE)
}


#Open files
data.counts <- read.table(
  opt$file,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

metadata <- read.delim(
  opt$metadata,
  header = TRUE,
  stringsAsFactors = F,
  sep = opt$delimiter,
  check.names = FALSE,
  row.names = 1
)

gene_metadata <- read.delim(
  opt$gene_metadata,
  header = TRUE,
  stringsAsFactors = F,
  sep = opt$gene_delimiter,
  check.names = FALSE,
  row.names = 1
)

#Retrieve significantly differentially expressed genes
DE_genes <- rownames(subset(gene_metadata, Significant == TRUE))

#Filter expression matrix
data <- data.counts[DE_genes,]

#Add signature score to expression matrix
data <- rbind(t(subset(metadata, select = c(opt$score))), data[,rownames(metadata)])

#Gene correlation
gene_corr <- rcorr(t(data), type = "pearson")

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
