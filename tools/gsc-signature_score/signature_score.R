#########################
#    Signature score    #
#########################

# Fourth step of the signature-based workflow
# Compute the signature score based on the
# geometric mean of the target gene expression
# and separate cells thanks to this signature
#score into 2 groups (high/low).

#Example of command
# Rscript 4-signature_score.R -f ../3-filter_genes/filterGenes.tsv -m ../1-filter_cells/filterCellsMetadata.tsv -l ARNT2,SALL2,SOX9,OLIG2,POU3F2 -o .

# Load necessary packages (install them if it's not the case)
requiredPackages = c('optparse', 'psych', 'ggplot2', 'gridExtra')
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
    c("-l", "--list"),
    default = NA,
    type = 'character',
    help = "List of genes comma separated"
  ),
  make_option(
    c("-p", "--percent"),
    default = 20,
    type = 'integer',
    help = "Percentage of dectection threshold [default : '%default' ]"
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

if (opt$file == "" | opt$list == "" & !(opt$help)) {
  stop("At least tw arguments must be supplied (count data --file option and gene list --list option).\n",
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

#Get vector of target genes
genes <- unlist(strsplit(opt$list, ","))

if (unique(genes %in% rownames(data.counts)) == F)
  stop("None of these genes are in your dataset: ", opt$list)

logical_genes <- rownames(data.counts) %in% genes

#Retrieve target genes in counts data
signature.counts <- subset(data.counts, logical_genes)


##Descriptive Statistics Function
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

#Remove poorly expressed genes from the signature 
kept_genes <- signature_stats$Percentage_Detection >= opt$percent

signature.counts <- signature.counts[kept_genes,]

#Add warnings
if (unique(kept_genes) == F) {
  stop(
    "None of these genes are detected in ",
    opt$percent,
    "% of your cells: ",
    rownames(signature_stats),
    ". You can be less stringent thanks to --percent parameter."
  )
}
if (length(unique(kept_genes)) > 1) {
  cat(
    "WARNINGS ! Following genes were removed from further analysis due to low gene expression :",
    paste(rownames(signature.counts)[!kept_genes], collapse = ",")
  )
}

#Replace 0 by 1 counts
signature.counts[signature.counts == 0] <- 1

#Geometric mean by cell
score <- apply(signature.counts, 2, geometric.mean)

#Add result in metadata
metadata_filtered <- metadata[names(score),]
metadata_filtered$Signature_category <- ifelse(score > mean(score), "HIGH", "LOW")
metadata_filtered$Signature_score <- score

#Re-arrange score matrix for plots
score <- data.frame(score = score,
                    order = rank(score, ties.method = "first"),
                    signature = metadata_filtered$Signature_category,
                    stringsAsFactors = F)

pdf(file = paste(opt$out, "signatureScore.pdf", sep = "/"))

ggplot(score, aes(x = order, y = score)) +
  geom_line() + 
  geom_segment(x = 0, xend = max(score$order[score$signature == "LOW"]), y = mean(score$score), yend = mean(score$score)) +
  geom_area(aes(fill = signature), alpha = .7) +
  scale_fill_manual(values=c("#ff0000", "#08661e")) +
  geom_text(aes(x = 1, y = mean(score)), label = "Mean", vjust = -0.3, colour = "black") +
  labs(title = "Score Signature", x = "Cell index", y = "Score")

density_score <- density(score$score)
ggplot(data.frame(density_score[1:2]), aes(x, y, fill = ifelse(x < mean(score$score), "LOW", "HIGH"))) +
  geom_line() +
  geom_vline(xintercept = mean(score$score)) +
  geom_text(x = mean(score$score), y = max(density_score$y), label = "Mean", hjust = -0.3, colour = "black") +
  geom_area(alpha = .7) +
  scale_fill_manual(values=c("#ff0000", "#08661e")) +
  ylim(0, max(density_score$y)) +
  labs(
    title = "Score density",
    x = paste("N =", density_score$n, "Bandwidth =", density_score$bw),
    y = "Density",
    fill = "Signature"
  )

#Check patient distribution in two groups
counts <- data.frame(table(metadata_filtered$Signature_category, metadata_filtered$Sample.name))

ggplot(data = counts, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = .8) +
  scale_fill_manual(values=c("#ff0000", "#08661e")) +
  labs(title = "Cell score distribution by patient", fill = "Score", x = "Patient")

#Check score independant of low expression
p_gene <- ggplot(metadata_filtered, aes(Signature_category, nGene)) +
  geom_violin(aes(fill = Signature_category), alpha = .5, trim = F, show.legend = F) +
  scale_fill_manual(values=c("#ff0000", "#08661e")) +
  geom_jitter() + labs(y = "Number of detected genes", x = "Signature")

p_counts <- ggplot(metadata_filtered, aes(Signature_category, total_counts)) +
  geom_violin(aes(fill = Signature_category), alpha = .5, trim = F, show.legend = F) +
  scale_fill_manual(values=c("#ff0000", "#08661e")) +
  geom_jitter() + labs(y = "Total counts", x = "Signature")

grid.arrange(p_gene, p_counts, ncol = 2, top = "Comparison of expression level based on their Signature score (HIGH vs LOW)")

dev.off()

#Save file
write.table(
  metadata_filtered,
  paste(opt$out, "signatureScoreMetadata.tsv", sep = "/"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)