# Performs multi-correlation analysis between the vectors of gene expressions
# in single cell RNAseq libraries and the vectors of signature scores in these
# same single cell RNAseq libraries.
# Example of command
# Rscript correlations_with_signature.R --expression_file <expression_data.tsv>
#                                       --signatures_file <signature_scores.tsv>
#                                       --sep "\t"
#                                       --colnames "T"
#                                       --gene_corr <gene-gene corr file>
#                                       --gene_corr_pval <gene-gene corr pvalues file>
#                                       --sig_corr <genes correlation to signature file>

options(show.error.messages = FALSE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library(optparse)
library(Hmisc)

# Arguments
option_list <- list(
  make_option(
    "--sep",
    default = "\t",
    type = "character",
    help = "File separator, must be the same for all input files [default : '%default' ]"
  ),
  make_option(
    "--colnames",
    default = TRUE,
    type = "logical",
    help = "Consider first lines as header (must stand for all input files) [default : '%default' ]"
  ),  
  make_option(
    "--expression_file",
    default = NA,
    type = "character",
    help = "Input file that contains log2(CPM +1) expression values"
  ),
  make_option(
    "--signatures_file",
    default = NA,
    type = "character",
    help = "Input file that contains cell signature"
  ),
  make_option(
    "--sig_corr",
    default = "sig_corr.tsv",
    type = "character",
    help = "signature correlations output [default : '%default' ]"
  ),
  make_option(
    "--gene_corr",
    default = "gene_corr.tsv",
    type = "character",
    help = "genes-genes correlations output [default : '%default' ]"
  ),
  make_option(
    "--gene_corr_pval",
    default = "gene_corr_pval.tsv",
    type = "character",
    help = "genes-genes correlations pvalues output [default : '%default' ]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

# Open files
data <- read.delim(
  opt$expression_file,
  header = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = FALSE
)
signature <- read.delim(
  opt$signatures_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1,
  sep = opt$sep,
  check.names = FALSE
)


# keep only signatures that are in the expression dataframe
signature <- subset(signature, rownames(signature) %in% colnames(data))

# Add signature score to expression matrix
data <- rbind(t(signature), data)

# Gene correlation
gene_corr <- rcorr(t(data), type = "pearson") # transpose because we correlate genes, not cells

# Gene correlation with signature score
gene_signature_corr <- cbind.data.frame(gene = colnames(gene_corr$r),
                                        Pearson_correlation = gene_corr$r[, 1], 
                                        p_value = gene_corr$P[, 1])
gene_signature_corr <- gene_signature_corr[ order(gene_signature_corr[, 2], decreasing = TRUE), ]


# Save files
write.table(
  gene_signature_corr,
  file = opt$sig_corr,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

r_genes <- data.frame(gene = rownames(gene_corr$r), gene_corr$r) # add rownames as a variable for output
p_genes <- data.frame(gene = rownames(gene_corr$P), gene_corr$P) # add rownames as a variable for output
write.table(
  r_genes[-1, -2],
  file = opt$gene_corr,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
write.table(
  p_genes[-1, -2], 
  file = opt$gene_corr_pval,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
