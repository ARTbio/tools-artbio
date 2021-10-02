# load packages that are provided in the conda env
options(show.error.messages = F,
       error = function() {
           cat(geterrmessage(), file = stderr()); q("no", 1, F)
           }
        )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)
library(rjson)
library(grid)
library(gridExtra)
library(scales)
library(RColorBrewer)

# Arguments
option_list <- list(
  make_option(
    "--inputs",
    default = NA,
    type = "character",
    help = "json formatted dictionary of datasets and their paths"
  ),
  make_option(
    "--genome",
    default = NA,
    type = "character",
    help = "genome name in the BSgenome bioconductor package"
  ),
  make_option(
    "--levels",
    default = NA,
    type = "character",
    help = "path to the tab separated file describing the levels in function of datasets"
  ),
  make_option(
    "--cosmic_version",
    default = "v2",
    type = "character",
    help = "Version of the Cosmic Signature set to be used to express mutational patterns"
  ),
  make_option(
    "--signum",
    default = 2,
    type = "integer",
    help = "selects the N most significant signatures in samples to express mutational patterns"
  ),
  make_option(
    "--nrun",
    default = 2,
    type = "integer",
    help = "Number of runs to fit signatures"
  ),
  make_option(
    "--rank",
    default = 2,
    type = "integer",
    help = "number of ranks to display for parameter optimization"
  ),
    make_option(
    "--newsignum",
    default = 2,
    type = "integer",
    help = "Number of new signatures to be captured"
  ),
  make_option(
    "--output_spectrum",
    default = NA,
    type = "character",
    help = "path to output dataset"
  ),
  make_option(
    "--output_denovo",
    default = NA,
    type = "character",
    help = "path to output dataset"
  ),
  make_option(
    "--sigmatrix",
    default = NA,
    type = "character",
    help = "path to signature matrix"
  ),
  make_option(
    "--output_cosmic",
    default = NA,
    type = "character",
    help = "path to output dataset"
  ),
  make_option(
    "--sig_contrib_matrix",
    default = NA,
    type = "character",
    help = "path to signature contribution matrix"
  ),
  make_option(
    c("-r", "--rdata"),
    type = "character",
    default = NULL,
    help = "Path to RData output file")
)

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

################ Manage input data ####################
json_dict <- opt$inputs
parser <- newJSONParser()
parser$addData(json_dict)
fileslist <- parser$getObject()
vcf_paths <- attr(fileslist, "names")
element_identifiers <- unname(unlist(fileslist))
ref_genome <- opt$genome
vcf_table <- data.frame(element_identifier = as.character(element_identifiers), path = vcf_paths)

library(MutationalPatterns)
library(ref_genome, character.only = TRUE)
library(ggplot2)

# Load the VCF files into a GRangesList:
vcfs <- read_vcfs_as_granges(vcf_paths, element_identifiers, ref_genome)
library(plyr)
if (!is.na(opt$levels)[1]) {  # manage levels if there are
    levels_table  <- read.delim(opt$levels, header = FALSE,
                                col.names = c("element_identifier", "level"))
    } else {
    levels_table <- data.frame(element_identifier = vcf_table$element_identifier,
                               level = rep("nolabels", length(vcf_table$element_identifier)))
}
metadata_table <- join(vcf_table, levels_table, by = "element_identifier")
tissue <- as.vector(metadata_table$level)
detach(package:plyr)

##### This is done for any section ######
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
col_vector <- col_vector[c(-32, -34, -39)] # 67-color palette

###### Section 1 Mutation characteristics and spectrums #############
if (!is.na(opt$output_spectrum)[1]) {
    pdf(opt$output_spectrum, paper = "special", width = 11.69, height = 11.69)
    type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

    # mutation spectrum, total or by sample

    if (length(levels(factor(levels_table$level))) == 1) {
        p1 <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE)
        plot(p1)
    } else {
        p2 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE) # by levels
        p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE) # total
        grid.arrange(p2, p3, ncol = 2, widths = c(4, 2.3), heights = c(4, 1))
   }
    plot_96_profile(mut_mat, condensed = TRUE)
    dev.off()
}

###### Section 2: De novo mutational signature extraction using NMF #######
# opt$rank cannot be higher than the number of samples and
# likewise, opt$signum cannot be higher thant the number of samples
if (!is.na(opt$output_denovo)[1]) {

    if (opt$rank > length(element_identifiers)) {
        opt$rank <- length(element_identifiers)
        }
    if (opt$signum > length(element_identifiers)) {
        opt$signum <- length(element_identifiers)
        }
    pseudo_mut_mat <- mut_mat + 0.0001 # First add a small pseudocount to the mutation count matrix
    # Use the NMF package to generate an estimate rank plot
    library("NMF")
    estimate <- nmf(pseudo_mut_mat, rank = 1:opt$rank, method = "brunet", nrun = opt$nrun, seed = 123456)
    # And plot it
    pdf(opt$output_denovo, paper = "special", width = 11.69, height = 11.69)
    p4 <- plot(estimate)
    grid.arrange(p4)
    # Extract 4 (PARAMETIZE) mutational signatures from the mutation count matrix with extract_signatures
    # (For larger datasets it is wise to perform more iterations by changing the nrun parameter
    # to achieve stability and avoid local minima)
    nmf_res <- extract_signatures(pseudo_mut_mat, rank = opt$newsignum, nrun = opt$nrun)
    # Assign signature names
    colnames(nmf_res$signatures) <- paste0("NewSig_", 1:opt$newsignum)
    rownames(nmf_res$contribution) <- paste0("NewSig_", 1:opt$newsignum)
    # Plot the 96-profile of the signatures:
    p5 <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
    new_sig_matrix <- reshape2::dcast(p5$data, substitution + context ~ sample, value.var = "freq")
    new_sig_matrix <- format(new_sig_matrix, scientific = TRUE)
    write.table(new_sig_matrix, file = opt$sigmatrix, quote = FALSE, row.names = FALSE, sep = "\t")
    grid.arrange(p5)
    # Visualize the contribution of the signatures in a barplot
    pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative", coord_flip = TRUE)
    # Visualize the contribution of the signatures in absolute number of mutations
    pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", coord_flip = TRUE)
    # Combine the two plots:
    grid.arrange(pc1, pc2)

    # The relative contribution of each signature for each sample can also be plotted as a heatmap with
    # plot_contribution_heatmap, which might be easier to interpret and compare than stacked barplots.
    # The samples can be hierarchically clustered based on their euclidean dis- tance. The signatures
    # can be plotted in a user-specified order.
    # Plot signature contribution as a heatmap with sample clustering dendrogram and a specified signature order:
    pch1 <- plot_contribution_heatmap(nmf_res$contribution,
                                      sig_order = paste0("NewSig_", 1:opt$newsignum))
    # Plot signature contribution as a heatmap without sample clustering:
    pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples = FALSE)
    #Combine the plots into one figure:
    grid.arrange(pch1, pch2, ncol = 2, widths = c(2, 1.6))

    # Compare the reconstructed mutational profile with the original mutational profile:
    plot_compare_profiles(pseudo_mut_mat[, 1],
                          nmf_res$reconstructed[, 1],
                          profile_names = c("Original", "Reconstructed"),
                          condensed = TRUE)
    dev.off()
    }
##### Section 3: Find optimal contribution of known signatures: COSMIC mutational signatures ####

if (!is.na(opt$output_cosmic)[1]) {
    pdf(opt$output_cosmic, paper = "special", width = 11.69, height = 11.69)
    pseudo_mut_mat <- mut_mat + 0.0001 # First add a small psuedocount to the mutation count matrix
    if (opt$cosmic_version == "v2") {
        sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
        cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
        new_order <- match(row.names(pseudo_mut_mat), cancer_signatures$Somatic.Mutation.Type)
        cancer_signatures <- cancer_signatures[as.vector(new_order), ]
        row.names(cancer_signatures) <- cancer_signatures$Somatic.Mutation.Type
        cancer_signatures <- as.matrix(cancer_signatures[, 4:33])
        colnames(cancer_signatures) <- gsub("Signature.", "", colnames(cancer_signatures)) # shorten signature labels
        cosmic_tag <- "Signatures (Cosmic v2, March 2015)"
        cosmic_colors <- col_vector[1:30]
        names(cosmic_colors) <- colnames(cancer_signatures)
        } else {
        sp_url <- "https://raw.githubusercontent.com/ARTbio/startbio/master/sigProfiler_SBS_signatures_2019_05_22.tsv"
        cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
        new_order <- match(row.names(pseudo_mut_mat), cancer_signatures$Somatic.Mutation.Type)
        cancer_signatures <- cancer_signatures[as.vector(new_order), ]
        row.names(cancer_signatures) <- cancer_signatures$Somatic.Mutation.Type
        cancer_signatures <- as.matrix(cancer_signatures[, 4:70])
        colnames(cancer_signatures) <- gsub("SBS", "", colnames(cancer_signatures)) # shorten signature labels
        cosmic_tag <- "Signatures (Cosmic v3, May 2019)"
        cosmic_colors <- col_vector[1:67]
        names(cosmic_colors) <- colnames(cancer_signatures)
        }

    # Plot mutational profiles of the COSMIC signatures
    if (opt$cosmic_version == "v2") {
        p6 <- plot_96_profile(cancer_signatures, condensed = TRUE, ymax = 0.3)
        grid.arrange(p6, top = textGrob("COSMIC signature profiles", gp = gpar(fontsize = 12, font = 3)))
    } else {
        p6 <- plot_96_profile(cancer_signatures[, 1:33], condensed = TRUE, ymax = 0.3)
        p6bis <- plot_96_profile(cancer_signatures[, 34:67], condensed = TRUE, ymax = 0.3)
        grid.arrange(p6, top = textGrob("COSMIC signature profiles (on two pages)",
                     gp = gpar(fontsize = 12, font = 3)))
        grid.arrange(p6bis, top = textGrob("COSMIC signature profiles (continued)",
                     gp = gpar(fontsize = 12, font = 3)))
    }


    # Find optimal contribution of COSMIC signatures to reconstruct 96 mutational profiles
    fit_res <- fit_to_signatures(pseudo_mut_mat, cancer_signatures)

    # Plot contribution barplots
    pc3 <- plot_contribution(fit_res$contribution, cancer_signatures, coord_flip = T, mode = "absolute")
    pc4 <- plot_contribution(fit_res$contribution, cancer_signatures, coord_flip = T, mode = "relative")
    if (is.na(opt$levels)[1]) {  # if there are NO levels to display in graphs
        pc3_data <- pc3$data
        pc3 <- ggplot(pc3_data, aes(x = Sample, y = Contribution, fill = as.factor(Signature))) +
               geom_bar(stat = "identity", position = "stack") +
               coord_flip() +
               scale_fill_manual(name = "Cosmic\nSignatures", values = cosmic_colors) +
               labs(x = "Samples", y = "Absolute contribution") + theme_bw() +
               theme(panel.grid.minor.x = element_blank(),
                     panel.grid.major.x = element_blank(),
                     legend.position = "right",
                     text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, hjust = 1))
        pc4_data <- pc4$data
        pc4 <- ggplot(pc4_data, aes(x = Sample, y = Contribution, fill = as.factor(Signature))) +
               geom_bar(stat = "identity", position = "fill") +
               coord_flip() +
               scale_fill_manual(name = "Cosmic\nSignatures", values = cosmic_colors) +
               scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
               labs(x = "Samples", y = "Relative contribution") + theme_bw() +
               theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "right",
                     text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, hjust = 1))
    }
    #####
    # ggplot2 alternative
    if (!is.na(opt$levels)[1]) {  # if there are levels to display in graphs
        pc3_data <- pc3$data
        pc3_data <- merge(pc3_data, metadata_table[, c(1, 3)], by.x = "Sample", by.y = "element_identifier")
        pc3 <- ggplot(pc3_data, aes(x = Sample, y = Contribution, fill = as.factor(Signature))) +
               geom_bar(stat = "identity", position = "stack") +
               scale_fill_manual(name = "Cosmic\nSignatures", values = cosmic_colors) +
               labs(x = "Samples", y = "Absolute contribution") + theme_bw() +
               theme(panel.grid.minor.x = element_blank(),
                     panel.grid.major.x = element_blank(),
                     legend.position = "right",
                     text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, hjust = 1)) +
               facet_grid(~level, scales = "free_x", space = "free")
        pc4_data <- pc4$data
        pc4_data <- merge(pc4_data, metadata_table[, c(1, 3)], by.x = "Sample", by.y = "element_identifier")
        pc4 <- ggplot(pc4_data, aes(x = Sample, y = Contribution, fill = as.factor(Signature))) +
               geom_bar(stat = "identity", position = "fill") +
               scale_fill_manual(name = "Cosmic\nSignatures", values = cosmic_colors) +
               scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
               labs(x = "Samples", y = "Relative contribution") + theme_bw() +
               theme(panel.grid.minor.x = element_blank(),
                     panel.grid.major.x = element_blank(),
                     legend.position = "right",
                     text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, hjust = 1)) +
               facet_grid(~level, scales = "free_x", space = "free")
    }
    # Combine the two plots:
    grid.arrange(pc3, pc4,
                 top = textGrob("Absolute and Relative Contributions of Cosmic signatures to mutational patterns",
                 gp = gpar(fontsize = 12, font = 3)))

    #### pie charts of comic signatures contributions in samples ###
    library(reshape2)
    library(dplyr)
    if (length(levels(factor(levels_table$level))) < 2) {
        fit_res_contrib <- as.data.frame(fit_res$contribution)
        worklist <- cbind(signature = rownames(fit_res$contribution),
                          level = rep("nolabels", length(fit_res_contrib[, 1])),
                          fit_res_contrib,
                          sum = rowSums(fit_res_contrib))
        worklist <- worklist[order(worklist[, "sum"], decreasing = T), ]
        worklist <- worklist[1:opt$signum, ]
        worklist <- worklist[, -length(worklist[1, ])]
        worklist <- melt(worklist)
        worklist <- worklist[, c(1, 3, 4, 2)]
    } else {
        worklist <- list()
        for (i in levels(factor(levels_table$level))) {
             fit_res$contribution[, levels_table$element_identifier[levels_table$level == i]] -> worklist[[i]]
             sum <- rowSums(as.data.frame(worklist[[i]]))
             worklist[[i]] <- cbind(worklist[[i]], sum)
             worklist[[i]] <- worklist[[i]][order(worklist[[i]][, "sum"], decreasing = T), ]
             worklist[[i]] <- worklist[[i]][1:opt$signum, ]
             worklist[[i]] <- worklist[[i]][, -length(as.data.frame(worklist[[i]]))]
        }
        worklist <- as.data.frame(melt(worklist))
        worklist[, 2] <- paste0(worklist[, 4], " - ", worklist[, 2])
    }

    colnames(worklist) <- c("signature", "sample", "value", "level")
    worklist <- as.data.frame(worklist %>% group_by(sample) %>% mutate(value = value / sum(value) * 100))
    worklist$pos <- cumsum(worklist$value) - worklist$value / 2
    worklist$label <- factor(worklist$signature)
    worklist$signature <- factor(worklist$signature)
    p7 <-  ggplot(worklist, aes(x = "", y = value, group = signature, fill = signature)) +
              geom_bar(width = 1, stat = "identity") +
              geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "black", size = 3) +
              coord_polar("y", start = 0) + facet_wrap(.~sample) +
              labs(x = "", y = "Samples", fill = cosmic_tag) +
              scale_fill_manual(name = paste0(opt$signum, " most contributing\nsignatures\n(in each label/tissue)"),
                                values = cosmic_colors[levels(worklist$signature)]) +
              theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid  = element_blank())
    grid.arrange(p7)

    # Plot relative contribution of the cancer signatures in each sample as a heatmap with sample clustering
    if (length(vcf_paths) > 1) {
        p8 <- plot_contribution_heatmap(fit_res$contribution, cluster_samples = TRUE, method = "complete")
        grid.arrange(p8)
    }

    # export relative contribution matrix
    if (!is.na(opt$sig_contrib_matrix)) {
        output_table <- t(fit_res$contribution) / rowSums(t(fit_res$contribution))
        colnames(output_table) <- paste0("s", colnames(output_table))
        if (length(levels(factor(levels_table$level))) > 1) {
            output_table <- data.frame(sample = paste0(metadata_table[metadata_table$element_identifier == colnames(fit_res$contribution),
                                                                    3], "-", colnames(fit_res$contribution)),
                                       output_table)
            } else {
        output_table <- data.frame(sample = rownames(output_table), output_table)
        }
        write.table(output_table, file = opt$sig_contrib_matrix, sep = "\t", quote = F, row.names = F)
    }

    # calculate all pairwise cosine similarities
    cos_sim_ori_rec <- cos_sim_matrix(pseudo_mut_mat, fit_res$reconstructed)
    # extract cosine similarities per sample between original and reconstructed
    cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))

    # We can use ggplot to make a barplot of the cosine similarities between the original and
    # reconstructed mutational profile of each sample. This clearly shows how well each mutational
    # profile can be reconstructed with the COSMIC mutational signatures. Two identical profiles
    # have a cosine similarity of 1. The lower the cosine similarity between original and
    # reconstructed, the less well the original mutational profile can be reconstructed with
    # the COSMIC signatures. You could use, for example, cosine similarity of 0.95 as a cutoff.

    # Adjust data frame for plotting with gpplot
    colnames(cos_sim_ori_rec) <- "cos_sim"
    cos_sim_ori_rec$sample <- row.names(cos_sim_ori_rec)
    # Make barplot
    p9 <- ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = sample)) +
                      geom_bar(stat = "identity", fill = "skyblue4") +
                      coord_cartesian(ylim = c(0.8, 1)) +
                      # coord_flip(ylim=c(0.8,1)) +
                      ylab("Cosine similarity\n original VS reconstructed") +
                      xlab("") +
                      # Reverse order of the samples such that first is up
                      # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
                      theme_bw() +
                      theme(panel.grid.minor.y = element_blank(),
                      panel.grid.major.y = element_blank()) +
                      # Add cut.off line
                      geom_hline(aes(yintercept = .95))
    grid.arrange(p9, top = textGrob("Similarity between true and reconstructed profiles (with all Cosmic sig.)", gp = gpar(fontsize = 12, font = 3)))
    dev.off()
}


# Output RData file
if (!is.null(opt$rdata)) {
  save.image(file = opt$rdata)
}
