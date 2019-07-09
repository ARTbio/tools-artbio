#########################
#    Signature score    #
#########################

# Compute the signature score based on the geometric mean of the target gene expression
# and split cells  in 2 groups (high/low) using this signature score.

# Example of command
# Rscript 4-signature_score.R --input <input.tsv> --genes  ARNT2,SALL2,SOX9,OLIG2,POU3F2
#                             --output <output.tab> --pdf <output.pdf>

# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()

library(optparse)
library(psych)
library(ggplot2)
library(gridExtra)

# Arguments
option_list = list(
  make_option(
    "--input",
    default = NA,
    type = 'character',
    help = "Input file that contains log2(CPM +1) values"
  ),
  make_option(
    "--sep",
    default = '\t',
    type = 'character',
    help = "File separator [default : '%default' ]"
  ),
  make_option(
    "--colnames",
    default = TRUE,
    type = 'logical',
    help = "Consider first line as header ? [default : '%default' ]"
  ),  
  make_option(
    "--genes",
    default = NA,
    type = 'character',
    help = "List of genes comma separated"
  ),
  make_option(
    "--percentile_threshold",
    default = 20,
    type = 'integer',
    help = "detection threshold to keep a gene in signature set [default : '%default' ]"
  ),
  make_option(
    "--output",
    default = "./output.tab",
    type = 'character',
    help = "Output path [default : '%default' ]"
  ),
  make_option(
    "--stats",
    default = "./statistics.tab",
    type = 'character',
    help = "statistics path [default : '%default' ]"
  ),
  make_option(
    "--correlations",
    default = "./correlations.tab",
    type = 'character',
    help = "Correlations between signature genes  [default : '%default' ]"
  ),
  make_option(
    "--covariances",
    default = "./statistics.tab",
    type = 'character',
    help = "Covariances between signature genes [default : '%default' ]"
  ),
  make_option(
    "--pdf",
    default = "~/output.pdf",
    type = 'character',
    help = "pdf path [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

# Take input data
data.counts <- read.table(
  opt$input,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

# Get vector of target genes
genes <- unlist(strsplit(opt$genes, ","))

if (length(unique(genes %in% rownames(data.counts))) == 1) {
  if (unique(genes %in% rownames(data.counts)) == F)
    stop("None of these genes are in your dataset: ", opt$genes)
}
    
logical_genes <- rownames(data.counts) %in% genes

# Retrieve target genes in counts data
signature.counts <- subset(data.counts, logical_genes)

# compute covariance
signature.covariances <- as.data.frame(cov(t(signature.counts)))
signature.covariances <- cbind(gene=rownames(signature.covariances), signature.covariances)
write.table(signature.covariances, file=opt$covariances, quote=F, row.names=F, sep="\t")

# compute signature.correlations
signature.correlations <- as.data.frame(cov(t(signature.counts)))
signature.correlations <- cbind(gene=rownames(signature.correlations), signature.correlations)
write.table(signature.correlations, file=opt$correlations, quote=F, row.names=F, sep="\t")

## Descriptive Statistics Function
descriptive_stats = function(InputData) {
  SummaryData = data.frame(
    mean = rowMeans(InputData),
    SD = apply(InputData, 1, sd),
    Variance = apply(InputData, 1, var),
    Percentage_Detection = apply(InputData, 1, function(x, y = InputData) {
      (sum(x != 0) / ncol(y)) * 100
    })
  )
  return(SummaryData)
}

signature_stats <- descriptive_stats(signature.counts)

# Find poorly detected genes from the signature 
kept_genes <- signature_stats$Percentage_Detection >= opt$percentile_threshold

# Add warnings
if (length(unique(kept_genes)) > 1) {
  cat(
    "WARNINGS ! Following genes were removed from further analysis due to low gene expression :",
    paste(paste(rownames(signature.counts)[!kept_genes], round(signature_stats$Percentage_Detection[!kept_genes], 2), sep = " : "), collapse = ", "),
    "\n"
  )
} else {
  if (unique(kept_genes) == F) {
    stop(
      "None of these genes are detected in ",
      opt$percent,
      "% of your cells: ",
      paste(rownames(signature_stats), collapse = ", "),
      ". You can be less stringent thanks to --percent parameter."
    )
  }
}

# Remove genes poorly detected in the dataset
signature.counts <- signature.counts[kept_genes,]
    
# Replace 0 by 1 counts
signature.counts[signature.counts == 0] <- 1

# Geometric mean by cell
score <- apply(signature.counts, 2, geometric.mean) # geometric.mean requires psych

# Add results in signature_output
signature_output <- data.frame(
                         cell = names(score),
                         score = score,
                         rate = ifelse(score > mean(score), "HIGH", "LOW"),
                         nGenes = colSums(data.counts != 0),
                         total_counts = colSums(data.counts)
                         )

# statistics of input genes, signature genes first lines
statistics.counts <- rbind(subset(data.counts, logical_genes),
                    subset(data.counts, !logical_genes))  
statistics <- descriptive_stats(statistics.counts)
statistics <- cbind(gene=rownames(statistics), statistics)



# Re-arrange score matrix for plots
score <- data.frame(score = score,
                    order = rank(score, ties.method = "first"),
                    signature = signature_output$rate,
                    stringsAsFactors = F)

pdf(file = opt$pdf)
myplot <- ggplot(signature_output, aes(x=rate, y=score)) +
                 geom_violin(aes(fill = rate), alpha = .5, trim = F, show.legend = F, cex=0.5) +
                 geom_abline(slope=0, intercept=mean(score$score), lwd=.5, color="red") +
                 scale_fill_manual(values=c("#ff0000", "#08661e")) +
                 geom_jitter(size=0.2) + labs(y = "Score", x = "Rate") +
                 annotate("text", x = 0.55, y = mean(score$score), cex = 3, vjust=1.5,
                           color="black", label = mean(score$score), parse = TRUE) +
                 labs(title = "Violin plots of Cell signature scores")

print(myplot)
dev.off()

# Save file
write.table(
  signature_output,
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

write.table(
  statistics,
  opt$stats,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)