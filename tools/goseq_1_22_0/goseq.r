sink(stdout(), type = "message")
suppressWarnings(suppressMessages(library(goseq)))
suppressWarnings(suppressMessages(library(optparse)))

option_list <- list(
    make_option(c("-d", "--dge_file"), type="character", help="Path to file with differential gene expression result"),
    make_option(c("-w","--wallenius_tab"), type="character", help="Path to output file with P-values estimated using wallenius distribution."),
    make_option(c("-s","--sampling_tab"), type="character", default=FALSE, help="Path to output file with P-values estimated using wallenius distribution."),
    make_option(c("-n","--nobias_tab"), type="character", default=FALSE, help="Path to output file with P-values estimated using wallenius distribution and no correction for gene length bias."),
    make_option(c("-l","--length_bias_plot"), type="character", default=FALSE, help="Path to length-bias plot."),
    make_option(c("-sw","--sample_vs_wallenius_plot"), type="character", default=FALSE, help="Path to plot comparing sampling with wallenius p-values."),
    make_option(c("-padj", "--p_adj_column"), type="integer",help="Column that contains p. adjust values"),
    make_option(c("-c", "--cutoff"), type="double",dest="p_adj_cutoff",
                help="Genes with p.adjust below cutoff are considered not differentially expressed and serve as control genes"),
    make_option(c("-r", "--repcnt"), type="integer", default=100, help="Number of repeats for sampling"),
    make_option(c("-lf", "--length_file"), type="character", default="FALSE", help = "Path to tabular file mapping gene id to length"),
    make_option(c("-g", "--genome"), type="character", help = "Genome [used for looking up correct gene length]"),
    make_option(c("-i", "--gene_id"), type="character", help="Gene ID of gene column in DGE file"),
    make_option(c("-cat", "--use_genes_without_cat"), default=FALSE, type="logical", help="A boolean to indicate whether genes without a categorie should still be used. For example, a large number of gene may have no GO term annotated.  If thisoption is set to FALSE, those genes will be ignored in the calculation of p-values(default behaviour).  If this option is set to TRUE, then these genes will count towards  the  total  number  of  genes  outside  the  category  being  tested  (default behaviour prior to version 1.15.2)."
)
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
dge_file = args$dge_file
p_adj_column = args$p_adj_colum
p_adj_cutoff = args$p_adj_cutoff
length_file = args$length_file
genome = args$genome
gene_id = args$gene_id
wallenius_tab = args$wallenius_tab
sampling_tab = args$sampling_tab
nobias_tab = args$nobias_tab
length_bias_plot = args$length_bias_plot
sample_vs_wallenius_plot = args$sample_vs_wallenius_plot
repcnt = args$repcnt
use_genes_without_cat = args$use_genes_without_cat

# format DE genes into vector suitable for use with goseq
dge_table = read.delim(dge_file, header = TRUE, sep="\t", check.names = FALSE)
genes = as.integer(dge_table[,p_adj_column]<p_adj_cutoff)
names(genes) = dge_table[,1] # Assuming first row contains gene names

# Get gene lengths
if (length_file != "FALSE" ) {
  length_table = read.delim(length_file, header=TRUE, sep="\t", check.names=FALSE)
  row.names(length_table) = length_table[,1]
  gene_lengths = length_table[names(genes),]$length
  } else {
  gene_lengths = getlength(names(genes), genome, gene_id)
  }

# Estimate PWF

pdf(length_bias_plot)
pwf=nullp(genes, genome, gene_id, gene_lengths)
message = dev.off()

# Fetch GO annotations:
go_map=getgo(names(genes), genome, gene_id, fetch.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"))

# wallenius approximation of p-values
GO.wall=goseq(pwf, genome, gene_id, use_genes_without_cat = use_genes_without_cat, gene2cat=go_map)

GO.nobias=goseq(pwf, genome, gene_id, method="Hypergeometric", use_genes_without_cat = use_genes_without_cat, gene2cat=go_map)

# Sampling distribution
if (repcnt > 0) {
  GO.samp=goseq(pwf,genome, gene_id, method="Sampling", repcnt=repcnt, use_genes_without_cat = use_genes_without_cat, gene2cat=go_map)
  # Compare sampling with wallenius
  pdf(sample_vs_wallenius_plot)
  plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
     abline(0,1,col=3,lty=2)
  message = dev.off()
  write.table(GO.samp, sampling_tab, sep="\t", row.names = FALSE, quote = FALSE)
}


write.table(GO.wall, wallenius_tab, sep="\t", row.names = FALSE, quote = FALSE)
write.table(GO.nobias, nobias_tab, sep="\t", row.names = FALSE, quote = FALSE)

sessionInfo()

# Use the following to get a list of supported genomes / gene ids

# write.table(supportedGenomes(), "available_genomes.tab", row.names = FALSE, quote=FALSE)
# write.table(supportedGeneIDs(), "supported_gene_ids.tab", row.name = FALSE, quote = FALSE)
# write.table(table.summary, "input_gene_count_matrix.tab", row.names = FALSE, quote = FALSE)
