####################
#   Differential   #
#     analysis     #
####################

# Fifth step of the signature-based workflow
# Perform a differential analysis between 2
# groups high/low.

# Example of command
# Rscript 5-differential_analysis.R -f ../3-filter_genes/filterGenes.tsv -m ../4-signature_score/signatureScoreMetadata.tsv --name Signature_category -o .

# Load necessary packages (install them if it's not the case)
requiredPackages = c('optparse', 'Hmisc')
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
    c("-n", "--name"),
    default = "signature",
    type = 'character',
    help = "Column name of signature category (HIGH/LOW) in metadata file [default : '%default' ]"
  ), 
  make_option(
    "--fdr",
    default = 0.01,
    type = 'numeric',
    help = "FDR threshold [default : '%default' ]"
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

#Create a logical named vector of whether or not the cell is signature "high"
high_cells <- setNames(metadata[,opt$name] == "HIGH", rownames(metadata))

## Mann-Whitney test (Two-sample Wilcoxon test)
MW_test <- data.frame(t(apply(data.counts, 1, function(x) {
  do.call("cbind", wilcox.test(x[high_cells], x[!high_cells]))[, 1:4]
})), stringsAsFactors = F)

MW_test[,1:3] <- apply(MW_test[,1:3], 2, as.numeric) #Change type "chr" to "num"

#Benjamini-Hochberg correction and significativity
MW_test$p.adjust <- p.adjust(MW_test$p.value, method = "BH" , n = nrow(MW_test))
MW_test$Critical.value <- (rank(MW_test$p.value) / nrow(MW_test)) * opt$fdr
MW_test$Significant <- MW_test$p.adjust < opt$fdr

##Descriptive Statistics Function
descriptive_stats = function(InputData) {
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
rownames(results) <- results$Row.names

#Save files
write.table(
  results[,-1],
  paste(opt$out, "diffAnalysisGeneMetadata.tsv", sep = "/"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)
