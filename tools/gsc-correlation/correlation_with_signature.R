# Performs multi-correlation analysis between the vectors of gene expressions
# in single cell RNAseq libraries and the vectors of signature scores in these
# same single cell RNAseq libraries.

# Example of command
# Rscript correlations_with_signature.R --expression_file <expression_data.tsv>
#                                       --signatures_file <signature_scores.tsv>
#                                       --DE_genes <differentially_express_genes.tsv>
#                                       --score_header "score" --sig_corr <sig corr file>
#                                       --gene_corr <gene-gene corr file>
#                                       --gene_corr_pval <gene-gene corr pvalues file>

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

# Arguments
option_list = list(
  make_option(
    "--sep",
    default = '\t',
    type = 'character',
    help = "File separator, must be the same for all input files [default : '%default' ]"
  ),
  make_option(
    "--colnames",
    default = TRUE,
    type = 'logical',
    help = "Consider first lines as header (must stand for all input files) [default : '%default' ]"
  ),  
  make_option(
    "--expression_file",
    default = NA,
    type = 'character',
    help = "Input file that contains log2(CPM +1) expression values"
  ),
  make_option(
    "--signatures_file",
    default = NA,
    type = 'character',
    help = "Input file that contains cell signature"
  ),
  make_option(
    "--score_header",
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
    "--sig_corr",
    default = "sig_corr.tsv",
    type = 'character',
    help = "signature correlations output [default : '%default' ]"
  ),
  make_option(
    "--gene_corr",
    default = "gene_corr.tsv",
    type = 'character',
    help = "genes-genes correlations output [default : '%default' ]"
  ),
  make_option(
    "--gene_corr_pval",
    default = "gene_corr_pval.tsv",
    type = 'character',
    help = "genes-genes correlations pvalues output [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

# Open files
data.counts <- read.table(
  opt$expression_file,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)
signature <- read.delim(
  opt$signatures_file,
  header = T,
  stringsAsFactors = F,
  sep = opt$sep,
  check.names = F,
  row.names = 1
)
diff_expressed_genes <- read.delim(
  opt$DE_genes,
  header = T,
  stringsAsFactors = F,
  sep = opt$sep,
  check.names = F,
  row.names = 1
)

# Retrieve significantly differentially expressed genes
# here a header name is imposed, to do :  select on padj !
DE_genes <- rownames(subset(diff_expressed_genes, Significant == TRUE))

# Filter expression matrix
data <- data.counts[DE_genes,]

# Add signature score to expression matrix
data <- rbind(t(subset(signature, select = c(opt$score_header))), data[,rownames(signature)])

# Gene correlation
gene_corr <- rcorr(t(data), type = "pearson") # transpose because we correlate genes, not cells

# Gene correlation with signature score
gene_corr_score <- cbind.data.frame(gene = colnames(gene_corr$r),
                                    r = gene_corr$r[, opt$score_header], 
                                    P = gene_corr$P[, opt$score_header])

# Heatmap
 heatmaply(
   subset(gene_corr_score, select = r),
   xlab = "Signature Score", 
   ylab = "Differentially Expressed Genes",
   main = "Correlation between Signature and significantly DE genes",
   dendrogram = "row",
   colors = RdBu,
   showticklabels = c(F,F),
   limits = c(-1,1),
   margins = c(40,NA,NA,NA),
   file = "./plot.html"
 )


# Save files
write.table(
  gene_corr_score,
  file = opt$sig_corr,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

r_genes <- data.frame(gene=rownames(gene_corr$r), gene_corr$r) # add rownames as a variable for output
p_genes <- data.frame(gene=rownames(gene_corr$P), gene_corr$P) # add rownames as a variable for output
write.table(
  r_genes,
  file = opt$gene_corr,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)
write.table(
  p_genes, 
  file = opt$gene_corr_pval,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)
