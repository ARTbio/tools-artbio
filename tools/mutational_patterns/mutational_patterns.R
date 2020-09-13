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
vcf_table <- data.frame(path=vcf_files, sample_name=sample_names)
metadata_table <- merge(vcf_table, levels_table, by.x=2, by.y=1)
levels <- metadata_table$level
muts = mutations_from_vcf(vcfs[[1]])
types = mut_type(vcfs[[1]])
context = mut_context(vcfs[[1]], ref_genome)
type_context = type_context(vcfs[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
# p1 <- plot_spectrum(type_occurrences)
# p2 <- plot_spectrum(type_occurrences, CT = TRUE)
# p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)
# 
# plot(p2)
# p4 <- plot_spectrum(type_occurrences, by = levels, CT = TRUE, legend = TRUE)
# palette <- c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
# p5 <- plot_spectrum(type_occurrences, CT=TRUE, legend=TRUE, colors=palette)
# 
# plot(p4)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
# plot_96_profile(mut_mat[,1:length(as.data.frame(mut_mat))], condensed = TRUE)
mut_mat <- mut_mat + 0.0001
# library("NMF")
# estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=100, seed=123456)
# plot(estimate)
# nmf_res <- extract_signatures(mut_mat, rank = 4, nrun = 100)
# colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
# rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
# plot_96_profile(nmf_res$signatures, condensed = TRUE)
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])
# plot_96_profile(cancer_signatures, condensed = TRUE, ymax = 0.3)
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
# plot(hclust_cosmic)
cos_sim(mut_mat[,1], cancer_signatures[,1])
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = TRUE)
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
threshold <- tail(sort(unlist(rowSums(fit_res$contribution), use.names = FALSE)), opt$signum)[1]
select <- which(rowSums(fit_res$contribution) >= threshold) # ensure opt$signum best signatures in samples are retained, the others discarded
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = T, mode = "absolute")
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = T, mode = "relative")
plot_contribution_heatmap(fit_res$contribution, cluster_samples = TRUE, method = "complete")


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
