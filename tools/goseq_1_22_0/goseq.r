sink(stdout(), type = "message")
library(goseq)
library(optparse)

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
    make_option(c("-g", "--genome"), type="character", help = "Genome [used for looking up correct gene length]"),
    make_option(c("-i", "--gene_id"), type="character", help="Gene ID of gene column in DGE file")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
dge_file = args$dge_file
p_adj_column = args$p_adj_colum
p_adj_cutoff = args$p_adj_cutoff
genome = args$genome
gene_id = args$gene_id
wallenius_tab = args$wallenius_tab
sampling_tab = args$sampling_tab
nobias_tab = args$nobias_tab
length_bias_plot = args$length_bias_plot
sample_vs_wallenius_plot = args$sample_vs_wallenius_plot
repcnt = args$repcnt


# format DE genes into vector suitable for use with goseq
dge_table = read.delim(dge_file, header = TRUE, sep="\t", check.names = FALSE)
genes = as.integer(dge_table[,p_adj_column]<p_adj_cutoff)
names(genes) = dge_table[,1] # Assuming first row contains gene names

# Estimate PWF

pdf(length_bias_plot)
pwf=nullp(genes, genome , gene_id)
dev.off()
# Null dstribution wallenius
GO.wall=goseq(pwf, genome, gene_id)

GO.nobias=goseq(pwf, genome, gene_id, method="Hypergeometric")

# Sampling dsitribution
GO.samp=goseq(pwf,genome, gene_id, method="Sampling",repcnt=repcnt)

# Compare sampling with wallenius
pdf(sample_vs_wallenius_plot)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)
dev.off()


write.table(GO.wall, wallenius_tab, sep="\t", row.names = FALSE, quote = FALSE)
write.table(GO.samp, sampling_tab, sep="\t", row.names = FALSE, quote = FALSE)
write.table(GO.nobias, nobias_tab, sep="\t", row.names = FALSE, quote = FALSE)

# Use the following to get a list of supported genomes / gene ids

# write.table(supportedGenomes(), "available_genomes.tab", row.names = FALSE, quote=FALSE)
# write.table(supportedGeneIDs(), "supported_gene_ids.tab", row.name = FALSE, quote = FALSE)
# write.table(table.summary, "input_gene_count_matrix.tab", row.names = FALSE, quote = FALSE)
