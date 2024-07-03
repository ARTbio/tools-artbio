if (length(commandArgs(TRUE)) == 0) {
    system("Rscript cpm_tpm_rpk.R -h", intern = FALSE)
    q("no")
}


# load packages that are provided in the conda env
options(show.error.messages = FALSE,
    error = function() {
        cat(geterrmessage(), file = stderr())
        q("no", 1, FALSE)
    }
)
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8") # nolint⁠
warnings()
library(optparse)
library(ggplot2)
library(reshape2)
library(Rtsne) # nolint⁠
library(ggfortify)



#Arguments
option_list <- list(
    make_option(
        c("-d", "--data"),
        default = NA,
        type = "character",
        help = "Input file that contains count values to transform"
    ),
    make_option(
        c("-t", "--type"),
        default = "cpm",
        type = "character",
        help = "Transformation type, either cpm, tpm, rpkm or none[default : '%default' ]"
    ),
    make_option(
        c("-s", "--sep"),
        default = "\t",
        type = "character",
        help = "File separator [default : '%default' ]"
    ),
    make_option(
        c("-c", "--colnames"),
        default = TRUE,
        type = "logical",
        help = "Consider first line as header ? [default : '%default' ]"
    ),
    make_option(
        c("-f", "--gene"),
        default = NA,
        type = "character",
        help = "Two column of gene length file"
    ),
    make_option(
        c("-a", "--gene_sep"),
        default = "\t",
        type = "character",
        help = "Gene length file separator [default : '%default' ]"
    ),
    make_option(
        c("-b", "--gene_header"),
        default = TRUE,
        type = "logical",
        help = "Consider first line of gene length as header ? [default : '%default' ]"
    ),
    make_option(
        c("-l", "--log"),
        default = FALSE,
        type = "logical",
        help = "Should be log transformed as well ? (log2(data +1)) [default : '%default' ]"
    ),
    make_option(
        c("-o", "--out"),
        default = "res.tab",
        type = "character",
        help = "Output name [default : '%default' ]"
    ),
    make_option(
        "--visu",
        default = FALSE,
        type = "logical",
        help = "performs T-SNE [default : '%default' ]"
    ),
    make_option(
        "--tsne_labels",
        default = FALSE,
        type = "logical",
        help = "add labels to t-SNE plot [default : '%default' ]"
    ),
    make_option(
        "--seed",
        default = 42,
        type = "integer",
        help = "Seed value for reproducibility [default : '%default' ]"
    ),
    make_option(
        "--perp",
        default = 5.0,
        type = "numeric",
        help = "perplexity [default : '%default' ]"
    ),
    make_option(
        "--theta",
        default = 1.0,
        type = "numeric",
        help = "theta [default : '%default' ]"
    ),
    make_option(
        c("-D", "--tsne_out"),
        default = "tsne.pdf",
        type = "character",
        help = "T-SNE pdf [default : '%default' ]"
    ),
    make_option(
        "--pca_out",
        default = "pca.pdf",
        type = "character",
        help = "PCA pdf [default : '%default' ]"
    )

)

opt <- parse_args(OptionParser(option_list = option_list),
                                    args = commandArgs(trailingOnly = TRUE))

if (opt$data == "" && !(opt$help)) {
    stop("At least one argument must be supplied (count data).\n",
             call. = FALSE)
} else if ((opt$type == "tpm" || opt$type == "rpkm") && opt$gene == "") {
    stop("At least two arguments must be supplied (count data and gene length file).\n",
             call. = FALSE)
} else if (opt$type != "tpm" && opt$type != "rpkm" && opt$type != "cpm" && opt$type != "none") {
    stop("Wrong transformation requested (--type option) must be : cpm, tpm or rpkm.\n",
             call. = FALSE)
}

if (opt$sep == "tab") {
    opt$sep <- "\t"
}
if (opt$gene_sep == "tab") {
    opt$gene_sep <- "\t"
}

cpm <- function(count) {
    (count / colSums(count)) * 1000000
}


rpk <- function(count, length) {
    count / (length / 1000)
}

tpm <- function(count, length) {
    rpk <- rpk(count, length)
    per_million_factor <- colSums(rpk) / 1000000
    tpm <- rpk / per_million_factor
    return(tpm)
}

rpkm <- function(count, length) {
    rpk <- rpk(count, length)
    per_million_factor <- colSum(count) / 1000000
    rpkm <- rpk / per_million_factor
    return(rpkm)
}

#### running code ####

data <- read.delim(
    opt$data,
    check.names = FALSE,
    header = opt$colnames,
    row.names = 1,
    sep = opt$sep
)

if (opt$type == "tpm" || opt$type == "rpkm") {
    gene_length <- as.data.frame(
        read.delim(
            opt$gene,
            header = opt$gene_header,
            row.names = 1,
            sep = opt$gene_sep
        )
    )
    gene_length <- as.data.frame(gene_length[match(rownames(data), rownames(gene_length)), ], rownames(data))
}


if (opt$type == "cpm")
    res <- cpm(data)
if (opt$type == "tpm")
    res <- as.data.frame(apply(data, 2, tpm, length = gene_length), row.names = rownames(data))
if (opt$type == "rpkm")
    res <- as.data.frame(apply(data, 2, rpkm, length = gene_length), row.names = rownames(data))
if (opt$type == "none")
    res <- data
colnames(res) <- colnames(data)

if (opt$log == TRUE) {
    res <- log2(res + 1)
}

if (opt$visu == TRUE) {
    df <- res
    # filter and transpose df for tsne and pca
    df <- df[rowSums(df) != 0, ] # remove lines without information (with only 0 counts)
    tdf <- t(df)
    # make tsne and plot results
    set.seed(opt$seed) ## Sets seed for reproducibility
    tsne_out <- Rtsne(tdf, perplexity = opt$perp, theta = opt$theta)
    embedding <- as.data.frame(tsne_out$Y)
    embedding$Class <- as.factor(sub("Class_", "", rownames(tdf)))
    gg_legend <- theme(legend.position = "none")
    ggplot(embedding, aes(x = V1, y = V2)) +
        geom_point(size = 1, color = "red") +
        gg_legend +
        xlab("") +
        ylab("") +
        ggtitle("t-SNE") +
        if (opt$tsne_labels == TRUE) {
            geom_text(aes(label = Class), hjust = -0.2, vjust = -0.5, size = 2.5, color = "darkblue")
        }
    ggsave(file = opt$tsne_out, device = "pdf")
    # make PCA and plot result with ggfortify (autoplot)
    tdf_pca <- prcomp(tdf, center = TRUE, scale. = TRUE)
    if (opt$tsne_labels == TRUE) {
        autoplot(tdf_pca, shape = FALSE, label = TRUE, label.size = 2.5, label.vjust = 1.2,
                         label.hjust = 1.2,
                         colour = "darkblue") +
            geom_point(size = 1, color = "red") +
            xlab(paste("PC1", summary(tdf_pca)$importance[2, 1] * 100, "%")) +
            ylab(paste("PC2", summary(tdf_pca)$importance[2, 2] * 100, "%")) +
            ggtitle("PCA")
        ggsave(file = opt$pca_out, device = "pdf")
    } else {
        autoplot(tdf_pca, shape = TRUE, colour = "darkblue") +
            geom_point(size = 1, color = "red") +
            xlab(paste("PC1", summary(tdf_pca)$importance[2, 1] * 100, "%")) +
            ylab(paste("PC2", summary(tdf_pca)$importance[2, 2] * 100, "%")) +
            ggtitle("PCA")
        ggsave(file = opt$pca_out, device = "pdf")
    }
}

# at this stage, we select numeric columns and round theirs values to 8 decimals for cleaner output
is_num <- sapply(res, is.numeric)
res[is_num] <- lapply(res[is_num], round, 8)

write.table(
    cbind(Features = rownames(res), res),
    opt$out,
    col.names = opt$colnames,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)
