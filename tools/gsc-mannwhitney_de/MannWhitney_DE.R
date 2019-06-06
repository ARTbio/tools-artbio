####################
#   Differential   #
#     analysis     #
####################

# Fifth step of the signature-based workflow
# Perform a differential analysis between 2
# groups high/low.

# Example of command
# Rscript MannWhitney_DE.R --input <input.tsv> --sep <tab> --colnames <TRUE> --metadata <signature.tsv> --names <rate> --fdr <0.01> --output <diff_analysis.tsv>

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
    "--metadata",
    default = NA,
    type = 'character',
    help = "Input file that contains cells metadata"
  ),
  make_option(
    "--name",
    default = "signature",
    type = 'character',
    help = "Column name of rate category in metadata file. It must be a vector of two categories only : 'HIGH' and 'LOW'  [default : '%default' ]"
  ), 
  make_option(
    "--fdr",
    default = 0.01,
    type = 'numeric',
    help = "FDR threshold [default : '%default' ]"
  ),
  make_option(
    "--output",
    default = "~",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$input == "" | opt$metadata == "" & !(opt$help)) {
  stop("At least two arguments must be supplied (count data --input option and cell metadata --metadata option).\n",
       call. = FALSE)
}

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}
if (opt$sep == "space") {opt$sep = " "}

#Open files
data.counts <- read.table(
  opt$input,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

metadata <- read.delim(
  opt$metadata,
  header = TRUE,
  stringsAsFactors = F,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

metadata <- subset(metadata, rownames(metadata) %in% colnames(data.counts))

# Create a logical named vector of whether or not the cell is signature "high"
high_cells <- setNames(metadata[,opt$name] == "HIGH", rownames(metadata))

## Mann-Whitney test (Two-sample Wilcoxon test)
MW_test <- data.frame(t(apply(data.counts, 1, function(x) {
  do.call("cbind", wilcox.test(x[names(high_cells)[high_cells]], x[names(high_cells)[!high_cells]]))[, 1:2]
})), stringsAsFactors = F)

MW_test[,1:3] <- apply(MW_test, 2, as.numeric) #Change type "chr" to "num"

# Benjamini-Hochberg correction and significativity
MW_test$p.adjust <- p.adjust(MW_test$p.value, method = "BH" , n = nrow(MW_test))
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
    mean_LOW = rowMeans(InputData[,!high_cells]),
    mean_HIGH = rowMeans(InputData[, high_cells])
  )
  SummaryData$fold_change = SummaryData$mean_HIGH / SummaryData$mean_LOW
  return(SummaryData)
}

gene_stats <- descriptive_stats(data.counts)

results <- merge(gene_stats, MW_test, by = "row.names")
colnames(results)[1] <- "genes"

# Save files
write.table(
  results,
  opt$output,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)
