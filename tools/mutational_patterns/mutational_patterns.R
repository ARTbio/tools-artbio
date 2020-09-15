# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)
library(rjson)
library(MutationalPatterns)
library(ggplot2)

# Arguments
option_list = list(
  make_option(
    "--inputs",
    default = NA,
    type = 'character',
    help = "json formatted dictionary of datasets and their paths"
  ),
  make_option(
    "--genome",
    default = NA,
    type = 'character',
    help = "genome name in the BSgenome bioconductor package"
  ),
  make_option(
    "--levels",
    default = NA,
    type = 'character',
    help = "path to the tab separated file describing the levels in function of datasets"
  ),
  make_option(
    "--signum",
    default = 2,
    type = 'integer',
    help = "selects the N most significant signatures in samples to express mutational patterns"
  ),
  make_option(
    "--output",
    default = NA,
    type = 'character',
    help = "path to output dataset"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

json_dict <- opt$inputs
parser <- newJSONParser()
parser$addData(json_dict)
fileslist <- parser$getObject()
vcf_files <- attr(fileslist, "names")
sample_names <- unname(unlist(fileslist))
pdf(opt$output, paper = "special", width = 11.69, height = 11.69)
ref_genome <- opt$genome
library(ref_genome, character.only = TRUE)

# Load the VCF files into a GRangesList:
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
levels_table  <- read.delim(opt$levels, header=FALSE, col.names=c("sample_name","level"))
vcf_table <- data.frame(sample_name=sample_names, path=vcf_files)
library(dplyr)
metadata_table <- left_join(vcf_table, levels_table, by = "sample_name")
tissue <- as.character(metadata_table$level)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

# only useful if only 1 level
p1 <- plot_spectrum(type_occurrences, CT = TRUE)
plot(p1)
print(metadata_table$level)
# mutation spectrum, total or by sample
p2 <- plot_spectrum(type_occurrences, by = metadata_table$level, CT=TRUE, legend=TRUE) # by sample
p3 <- plot_spectrum(type_occurrences, CT=TRUE, legend=TRUE) # total
library("gridExtra")
grid.arrange(p2, p3, ncol=2, widths=c(4,2.3), heights=c(4,1))

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat, condensed = TRUE)


## De novo mutational signature extraction using NMF ##
mut_mat <- mut_mat + 0.0001 # First add a small psuedocount to the mutation count matrix
# Use the NMF package to generate an estimate rank plot
library("NMF")
estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)  # PARAMETERIZE nrun  and rank!
# And plot it
plot(estimate)
# Extract 4 (PARAMETIZE) mutational signatures from the mutation count matrix with extract_signatures
# (For larger datasets it is wise to perform more iterations by changing the nrun parameter
# to achieve stability and avoid local minima)
nmf_res <- extract_signatures(mut_mat, rank = 4, nrun = 20)  # PARAMETERIZE nrun  and rank!
# Assign signature names
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")  # PARAMETERIZE
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")  # PARAMETERIZE
# Plot the 96-profile of the signatures:
plot_96_profile(nmf_res$signatures, condensed = TRUE)
# Visualize the contribution of the signatures in a barplot
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode="relative", coord_flip = TRUE)
# Visualize the contribution of the signatures in absolute number of mutations
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode="absolute", coord_flip = TRUE)
# Combine the two plots:
grid.arrange(pc1, pc2)


# The relative contribution of each signature for each sample can also be plotted as a heatmap with
# plot_contribution_heatmap, which might be easier to interpret and compare than stacked barplots.
# The samples can be hierarchically clustered based on their euclidean dis- tance. The signatures
# can be plotted in a user-specified order.
# Plot signature contribution as a heatmap with sample clustering dendrogram and a specified signature order:
pch1 <- plot_contribution_heatmap(nmf_res$contribution,
                                  sig_order = c("Signature B", "Signature A", "Signature C", "Signature D"))  # PARAMETERIZE
# Plot signature contribution as a heatmap without sample clustering:
pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)
#Combine the plots into one figure:
grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6))

# Compare the reconstructed mutational profile with the original mutational profile:   
plot_compare_profiles(mut_mat[,1],
                      nmf_res$reconstructed[,1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE)

##### Find optimal contribution of known signatures: COSMIC mutational signatures ####


sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

# Plot mutational profile the COSMIC signatures
plot_96_profile(cancer_signatures, condensed = TRUE, ymax = 0.3)
# Hierarchically cluster the COSMIC signatures based on their similarity with average linkage
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)

# Similarity between mutational profiles and COSMIC signatures


# The similarity between each mutational profile and each COSMIC signature, can be calculated
# with cos_sim_matrix, and visualized with plot_cosine_heatmap. The cosine similarity reflects
# how well each mutational profile can be explained by each signature individually. The advantage
# of this heatmap representation is that it shows in a glance the similarity in mutational
# profiles between samples, while at the same time providing information on which signatures
# are most prominent. The samples can be hierarchically clustered in plot_cosine_heatmap.
# The cosine similarity between two mutational profiles/signatures can be calculated with cos_sim :
# cos_sim(mut_mat[,1], cancer_signatures[,1])

# Calculate pairwise cosine similarity between mutational profiles and COSMIC signatures
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
# Plot heatmap with specified signature order
plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = TRUE)

# Find optimal contribution of COSMIC signatures to reconstruct 96 mutational profiles

fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
# Select signatures with some contribution (above a threshold)
threshold <- tail(sort(unlist(rowSums(fit_res$contribution), use.names = FALSE)), opt$signum)[1]
select <- which(rowSums(fit_res$contribution) >= threshold) # ensure opt$signum best signatures in samples are retained, the others discarded
# Plot contribution barplots
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = T, mode = "absolute")
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = T, mode = "relative")
# Plot relative contribution of the cancer signatures in each sample as a heatmap with sample clustering
plot_contribution_heatmap(fit_res$contribution, cluster_samples = TRUE, method = "complete")
# Compare the reconstructed mutational profile of sample 1 with its original mutational profile
plot_compare_profiles(mut_mat[,1], fit_res$reconstructed[,1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE)

# Calculate the cosine similarity between all original and reconstructed mutational profiles with
# `cos_sim_matrix`

# calculate all pairwise cosine similarities
cos_sim_ori_rec <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
# extract cosine similarities per sample between original and reconstructed
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))

# We can use ggplot to make a barplot of the cosine similarities between the original and
# reconstructed mutational profile of each sample. This clearly shows how well each mutational
# profile can be reconstructed with the COSMIC mutational signatures. Two identical profiles
# have a cosine similarity of 1. The lower the cosine similarity between original and
# reconstructed, the less well the original mutational profile can be reconstructed with
# the COSMIC signatures. You could use, for example, cosine similarity of 0.95 as a cutoff.

# Adjust data frame for plotting with gpplot
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
# ggplot2 is already loaded
# Make barplot
ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
     geom_bar(stat="identity", fill = "skyblue4") +
     coord_cartesian(ylim=c(0.8, 1)) +
     # coord_flip(ylim=c(0.8,1)) +
     ylab("Cosine similarity\n original VS reconstructed") +
     xlab("") +
     # Reverse order of the samples such that first is up
     # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
     theme_bw() +
     theme(panel.grid.minor.y=element_blank(),
     panel.grid.major.y=element_blank()) +
     # Add cut.off line
     geom_hline(aes(yintercept=.95))


## pie charts for comic signatures in samples

sig5data <- as.data.frame(t(head(fit_res$contribution[select,])))
colnames(sig5data) <- gsub("nature", "", colnames(sig5data))
sig5data_percents <- sig5data / (apply(sig5data,1,sum)) * 100
sig5data_percents$sample <- rownames(sig5data_percents)
library(reshape2)
melted_sig5data_percents <-melt(data=sig5data_percents)
melted_sig5data_percents$label <- sub("Sig.", "", melted_sig5data_percents$variable)
melted_sig5data_percents$pos <- cumsum(melted_sig5data_percents$value) - melted_sig5data_percents$value/2
ggplot(melted_sig5data_percents, aes(x="", y=value, group=variable, fill=variable)) +
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color="black", size=3) +
  coord_polar("y", start=0) + facet_wrap(~ sample) +
  labs(x="", y="Samples", fill = "Signatures (Cosmic_v2,March 2015)") +
    theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
dev.off()
