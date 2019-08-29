####################
#   Differential   #
#     analysis     #
####################

# Perform a differential analysis between 2
# groups of cells.

# Example of command
# Rscript MannWhitney_DE.R --input <input.tsv> --sep <tab> --colnames <TRUE> --metadata <signature.tsv> --column_name <rate> --fdr <0.01> --output <diff_analysis.tsv>

# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)

#Arguments
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
    "--comparison_factor_file",
    default = NA,
    type = 'character',
    help = " A two column table : cell identifiers and a comparison factor that split cells in two categories (high/low, HOM/HET,...)"
  ),
  make_option(
    "--factor1",
    type = 'character',
    help = "level associated to the control condition in the factor file"
  ), 
  make_option(
    "--factor2",
    type = 'character',
    help = "level associated to the test condition in the factor file"
  ),
  make_option(
    "--fdr",
    default = 0.01,
    type = 'numeric',
    help = "FDR threshold [default : '%default' ]"
  ),
  make_option(
    "--log",
    default=FALSE,
    action="store_true",
    type = 'logical',
    help = "Expression data are log-transformed [default : '%default' ]"
  ),
  make_option(
    "--output",
    default = "results.tsv",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

#Open files
data.counts <- read.table(
  opt$input,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

metadata <- read.table(
  opt$comparison_factor_file,
  header = TRUE,
  stringsAsFactors = F,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

metadata <- subset(metadata, rownames(metadata) %in% colnames(data.counts))

# Create two logical named vectors for each factor level of cell signature
factor1_cells <- setNames(metadata[,1] == opt$factor1, rownames(metadata))
factor2_cells <- setNames(metadata[,1] == opt$factor2, rownames(metadata))

## Mann-Whitney test (Two-sample Wilcoxon test)
MW_test <- data.frame(t(apply(data.counts, 1, function(x) {
  do.call("cbind", wilcox.test(x[names(factor1_cells)[factor1_cells]], x[names(factor2_cells)[factor2_cells]]))[, 1:2]
})), stringsAsFactors = F)

# Benjamini-Hochberg correction and significativity
MW_test$p.adjust <- p.adjust(as.numeric(MW_test$p.value), method = "BH" , n = nrow(MW_test))
# MW_test$Critical.value <- (rank(MW_test$p.value) / nrow(MW_test)) * opt$fdr
MW_test$Significant <- MW_test$p.adjust < opt$fdr

## Descriptive Statistics Function
descriptive_stats <- function(InputData) {
  SummaryData = data.frame(
    mean = rowMeans(InputData),
    SD = apply(InputData, 1, sd),
    Variance = apply(InputData, 1, var),
    Percentage_Detection = apply(InputData, 1, function(x, y = InputData) {
      (sum(x != 0) / ncol(y)) * 100
    }),
    mean_condition2 = rowMeans(InputData[,factor2_cells]),
    mean_condition1 = rowMeans(InputData[, factor1_cells])
  )
  if(opt$log) {
  SummaryData$log2FC <- SummaryData$mean_condition2 - SummaryData$mean_condition1
  } else {
  SummaryData$log2FC <- log2(SummaryData$mean_condition2 / SummaryData$mean_condition1)
  }
  return(SummaryData)
}

gene_stats <- descriptive_stats(data.counts)

results <- merge(gene_stats, MW_test, by = "row.names")
colnames(results)[1] <- "genes"

## Annotate Significant column
results$Significant[results$Significant == T & !is.na(results$Significant)] <- ifelse(subset(results, Significant == T)$log2FC > 0, "UP", "DOWN")
results$Significant[results$Significant == F & !is.na(results$Significant)] <- "NS"


# Save files
write.table(
  results[order(results$p.adjust),],
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)
