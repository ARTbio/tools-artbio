##################################################################################################
# Running PATHIFIER (Drier et al., 2013)
# Based on the work of Author: Miguel Angel Garcia-Campos - Github: https://github.com/AngelCampos
##################################################################################################


options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr()); q("no", 1, F)
  }
)
# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library(pathifier)
library(optparse)
library(pheatmap)
library(scatterplot3d)
library(circlize)

option_list <- list(
  make_option(
    "--exp",
    type = "character",
    help = "Expression matrix"),
  make_option(
    "--sep",
    type = "character",
    default = "\t",
    help = "File separator [default : '%default']"
  ),
  make_option(
    "--genes",
    type = "character",
    help = "Gene sets Pathways : gmt format (one pathway per line : Name, description, genes (one by column), tab separated)"),
  make_option(
    "--is_normal",
    default = F,
    type = "logical",
    help = "Define normals cells in your data"),
  make_option(
    "--normals",
    type = "character",
    help = "A vector of sample status : 1 = Healthy, 0 = Tumor. Must be in the same order as in expression data"),
  make_option(
    "--logfile",
    type = "character",
    default = "log.txt",
    help = "Log file name [default : '%default']"
  ),
  make_option(
    "--max_stability",
    type = "logical",
    default = T,
    help = "If true, throw away components leading to low stability of sampling noise [default : '%default']"
  ),
  make_option(
    "--attempts",
    type = "integer",
    default = 10,
    help = "Number of runs to determine stability. [default : '%default']"
  ),
  make_option(
    "--min_std",
    type = "character",
    default = "0.4",
    help = "Minimum of standard deviation to filter out low variable genes.
    Use --min.std data, to use the minimum std of your data [default : '%default']"
  ),
  make_option(
    "--min_exp",
    type = "character",
    default = "4",
    help = "Minimum of gene expression to filter out low expressed genes.
    Use --min.exp data, to use the minimum expression of your data [default : '%default']"
  ),
  make_option(
    "--pds",
    type = "character",
    default = "PDS.tsv",
    help = "Output PDS (Pathway deregulation score) of Pathifier in tabular file [default : '%default']"
  ),
  make_option(
    "--heatmap_cluster_cells",
    type = "logical",
    default = TRUE,
    help = "Cluster columns (cells) in the heatmap [default : '%default']"
  ),
  make_option(
    "--heatmap_cluster_pathways",
    type = "logical",
    default = TRUE,
    help = "Cluster rows (pathways) in the heatmap [default : '%default']"
  ),
  make_option(
    "--heatmap_show_cell_labels",
    type = "logical",
    default = FALSE,
    help = "Print column names (cells) on the heatmap [default : '%default']"
  ),
  make_option(
    "--heatmap_show_pathway_labels",
    type = "logical",
    default = FALSE,
    help = "Print row names (pathways) on the heatmap [default : '%default']"
  ),
  make_option(
    "--plot",
    type = "character",
    default = "./plot.pdf",
    help = "Pathifier visualization [default : '%default']"
  ),
  make_option(
    "--rdata",
    type = "character",
    default = "./results.rdata",
    help = "Pathifier object (S4) [default : '%default']"))
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)
if (args$sep == "tab") {
  args$sep <- "\t"
}


# set seed for reproducibility
set.seed(123)

# Load expression data for PATHIFIER
exp_matrix <- as.matrix(read.delim(file = args$exp,
                                   as.is = T,
                                   row.names = 1,
                                   sep = args$sep,
                                   check.names = F))

# Load Genesets annotation
gene_sets_file <- file(args$genes, open = "r")
gene_sets <- readLines(gene_sets_file)
close(gene_sets_file)

# Generate a list that contains genes in genesets
gs <- strsplit(gene_sets, "\t")
names(gs) <- lapply(gs, function(x) x[1])
gs <- lapply(gs, function(x) x[-c(1:2)])

# Generate a list that contains the names of the genesets used
pathwaynames <- names(gs)

# Prepare data and parameters ##################################################
# Extract information from binary phenotypes. 1 = Normal, 0 = Tumor
if (args$is_normal == T) {
  normals <- read.delim(file = args$normals, h = F)
  normals <- as.logical(normals[, 2])
  n_exp_matrix <- exp_matrix[, normals]
} else n_exp_matrix <- exp_matrix

# Calculate MIN_STD
rsd <- apply(n_exp_matrix, 1, sd)
min_std <- quantile(rsd, 0.25)

# Calculate MIN_EXP
min_exp <- quantile(as.vector(as.matrix(exp_matrix)),
                    0.1) # Percentile 10 of data

# Filter low value genes. At least 10% of samples with values over min_exp
# Set expression levels < MIN_EXP to MIN_EXP
over <- apply(exp_matrix, 1, function(x) x > min_exp)
g_over <- apply(over, 2, mean)
g_over <- names(g_over)[g_over > 0.1]
exp_matrix_filtered <- exp_matrix[g_over, ]
exp_matrix_filtered[exp_matrix_filtered < min_exp] <- min_exp

# Set maximum 5000 genes with more variance
variable_genes <- names(sort(apply(exp_matrix_filtered, 1, var), decreasing = T))[1:5000]
variable_genes <- variable_genes[!is.na(variable_genes)]
exp_matrix_filtered <- exp_matrix_filtered[variable_genes, ]
allgenes <- as.vector(rownames(exp_matrix_filtered))


if (args$min_std == "data") {
  args$min_std <- min_std
} else args$min_std <- as.numeric(args$min_std)

if (args$min_exp == "data") {
  args$min_exp <- min_exp
} else args$min_exp <- as.numeric(args$min_exp)


# Open pdf
pdf(args$plot)

# Construct continuous color scale
col_score_fun <- colorRamp2(c(0, 0.5, 1), c("#4575B4", "#FFFFBF", "#D73027"))

# Run Pathifier
if (args$is_normal == T) {
  pds <- quantify_pathways_deregulation(exp_matrix_filtered,
                                        allgenes,
                                        gs,
                                        pathwaynames,
                                        normals,
                                        maximize_stability = args$max_stability,
                                        attempts = args$attempts,
                                        logfile = args$logfile,
                                        min_std = args$min_std,
                                        min_exp = args$min_exp)
  for (i in pathwaynames) {
    df <- data.frame(pds$curves[[i]][, 1:3],
                     normal = normals,
                     PDS = as.numeric(pds$scores[[i]]),
                     curve_order = as.vector(pds$curves_order[[i]]))
    ordered <- df[df$curve_order, ]


    layout(cbind(1:2, 1:2), heights = c(7, 1))
    sc3 <- scatterplot3d(ordered[, 1:3],
                         main = paste("Principal curve of", i),
                         box = F, pch = 19, type = "l")
    sc3$points3d(ordered[, 1:3], box = F, pch = 19,
                 col = col_score_fun(ordered$PDS))

    # Plot color scale legend
    par(mar = c(5, 3, 0, 3))
    plot(seq(min(ordered$PDS), max(ordered$PDS), length = 100), rep(0, 100), pch = 15,
        axes = TRUE, yaxt = "n", xlab = "Color scale of PDS", ylab = "", bty = "n",
        col = col_score_fun(seq(min(ordered$PDS), max(ordered$PDS), length = 100)))


    cols_status <- ifelse(ordered$normal, "blue", "red")
    sc3 <- scatterplot3d(ordered[, 1:3],
                         main = paste("Principal curve of", i),
                         box = F, pch = "", type = "l")
    sc3$points3d(ordered[, 1:3], box = F,
                 pch = ifelse(ordered$normal, 19, 8),
                 col = ifelse(ordered$normal, "blue", "red"))
    legend("topright", pch = c(19, 8), yjust = 0,
           legend = c("normal", "cancer"),
           col = c("blue", "red"), cex = 1.1)

    ## annotation for heatmap
    sample_status <- data.frame(Status = factor(ifelse(df$normal, "normal", "tumor")))
    rownames(sample_status) <- colnames(exp_matrix_filtered)
    color_status_heatmap <- list(Status = c(normal = "blue", tumor = "red"))
  }
} else{
  pds <- quantify_pathways_deregulation(exp_matrix_filtered,
                                        allgenes,
                                        gs,
                                        pathwaynames,
                                        maximize_stability = args$max_stability,
                                        attempts = args$attempts,
                                        logfile = args$logfile,
                                        min_std = args$min_std,
                                        min_exp = args$min_exp)
  for (i in pathwaynames) {
    df <- data.frame(pds$curves[[i]][, 1:3],
                     PDS = as.numeric(pds$scores[[i]]),
                     curve_order = as.vector(pds$curves_order[[i]]))
    ordered <- df[df$curve_order, ]

    layout(cbind(1:2, 1:2), heights = c(7, 1))
    sc3 <- scatterplot3d(ordered[, 1:3],
                         main = paste("Principal curve of", i),
                         box = F, pch = 19, type = "l")
    sc3$points3d(ordered[, 1:3], box = F, pch = 19,
                 col = col_score_fun(ordered$PDS))

    # Plot color scale legend
    par(mar = c(5, 3, 0, 3))
    plot(seq(min(ordered$PDS), max(ordered$PDS), length = 100), rep(0, 100), pch = 15,
        axes = TRUE, yaxt = "n", xlab = "Color scale of PDS", ylab = "", bty = "n",
        col = col_score_fun(seq(min(ordered$PDS), max(ordered$PDS), length = 100)))


  ## annotation for heatmap (for the moment none for this situation)
  sample_status <- NA
  color_status_heatmap <- NA
  }
}

## Create dataframe from Pathifier list and round score to 4 digits
pds_scores <- mapply(FUN = function(x) cbind(round(x, 4)), pds$scores)
dimnames(pds_scores) <- list(colnames(exp_matrix_filtered), names(pds$scores))

## plot heatmap
if (ncol(pds_scores) > 1) {
    pheatmap(t(pds_scores),
             main = "Heatmap of Pathway Deregulation Score",   # heat map title
             cluster_rows = args$heatmap_cluster_pathways,     # apply clustering method
             cluster_cols = args$heatmap_cluster_cells,        # apply clustering method

             #Additional Options
             ## Color labeled columns
             annotation_col = sample_status,
             annotation_colors = color_status_heatmap,
             show_rownames = args$heatmap_show_pathway_labels,
             show_colnames = args$heatmap_show_cell_labels,
             border_color = NA,
             legend = TRUE)
}
dev.off()


## write table
write.table(pds_scores,
            args$pds,
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

## write S4 pathifier object
save(pds, file = args$rdata)
