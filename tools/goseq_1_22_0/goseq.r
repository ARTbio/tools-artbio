biocLite("goseq")
library(edgeR)
library(goseq)

# Vars:
dge_file = "${dge_file}"
p_adj_column = "${p_adj_column}"
p_adj_cutoff = "${p_adj_cutoff}"
genome = "${genome}"
gene_id = "${gene_id}"
wallenius_tab = "${wallenius_tab}"
sampling_tab = "${sampling_tab}"
nobias_tab = "${nobias_tab}"
length_bias_plot = "${length_bias_plot}"
sample_vs_wallenius_plot = "${sample_vs_wallenius_plot}"
repcnt = "${random_sampling_repititions}"


# format DE genes into vector suitable for use with goseq
dge_table = read.delim(dge_file, header = FALSE, sep="\t", quote = FALSE, check.names = FALSE)
genes = as.integer(dge_table[,p_adj_column]<p_adj_cutoff)
names(genes) = row.names(dge_table[,1]) # Assuming first row contains gene names

# Estimate PWF

pdf(length_bias_plot)
pwf=nullp(genes, genome , gene_id)
dev.off()
# Null dstribution wallenius
GO.wall=goseq(pwf, genome, gene_id)

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

# Use the following to get a lsit of supported genomes / gene ids

# write.table(supportedGenomes(), "available_genomes.tab", row.names = FALSE, quote=FALSE)
# write.table(supportedGeneIDs(), "supported_gene_ids.tab", row.name = FALSE, quote = FALSE)
# write.table(table.summary, "input_gene_count_matrix.tab", row.names = FALSE, quote = FALSE)
