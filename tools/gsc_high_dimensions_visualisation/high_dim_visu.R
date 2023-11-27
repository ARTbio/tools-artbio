options(show.error.messages = FALSE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()

library(optparse)
library(FactoMineR)
library(factoextra)
library(Rtsne)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(ClusterR)
library(data.table)
library(Polychrome)

#### Arguments ####
option_list <- list(
  make_option(
    "--data",
    default = NA,
    type = "character",
    help = "Input file that contains expression value to visualise"
  ),
  make_option(
    "--labels",
    default = FALSE,
    type = "logical",
    help = "add labels in scatter plots [default : '%default' ]"
  ),
  make_option(
    "--factor",
    default = "",
    type = "character",
    help = "A two column table that specifies factor levels for contrasting data [default : '%default' ]"
  ),
  make_option(
    "--visu_choice",
    default = "PCA",
    type = "character",
    help = "visualisation method ('PCA', 'tSNE', 'HCPC') [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_seed",
    default = 42,
    type = "integer",
    help = "Seed value for reproducibility [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_dims",
    default = 2,
    type = "integer",
    help = "Output dimensionality [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_initial_dims",
    default = 50,
    type = "integer",
    help = "The number of dimensions that should be retained in the initial PCA step [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_perplexity",
    default = 5.0,
    type = "numeric",
    help = "perplexity [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_theta",
    default = 1.0,
    type = "numeric",
    help = "theta [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_max_iter",
    default = 1000,
    type = "integer",
    help = "max_iter [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_pca",
    default = TRUE,
    type = "logical",
    help = "Whether an initial PCA step should be performed [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_pca_center",
    default = TRUE,
    type = "logical",
    help = "Should data be centered before pca is applied? [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_pca_scale",
    default = FALSE,
    type = "logical",
    help = "Should data be scaled before pca is applied? [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_normalize",
    default = TRUE,
    type = "logical",
    help = "Should data be normalized internally prior to distance calculations? [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_exaggeration_factor",
    default = 12.0,
    type = "numeric",
    help = " Exaggeration factor used to multiply the P matrix in the first part of the optimization [default : '%default' ]"
  ),
  make_option(
    "--PCA_npc",
    default = 5,
    type = "integer",
    help = "number of dimensions kept in the results [default : '%default' ]"
  ),
  make_option(
    "--item_size",
    default = 1,
    type = "numeric",
    help = "Size of points/labels in PCA [default : '%default' ]"
  ),
  make_option(
    "--x_axis",
    default = 1,
    type = "integer",
    help = "PC to plot in the x axis [default : '%default' ]"
  ),
  make_option(
    "--y_axis",
    default = 2,
    type = "integer",
    help = "PC to plot in the y axis [default : '%default' ]"
  ),
  make_option(
    "--HCPC_ncluster",
    default = -1,
    type = "numeric",
    help = "nb.clust, number of clusters to consider in the hierarchical clustering. [default : -1 let HCPC to optimize the number]"
  ),
  make_option(
    "--HCPC_npc",
    default = 5,
    type = "integer",
    help = "npc, number of dimensions which are kept for HCPC analysis [default : '%default' ]"
  ),
  make_option(
    "--HCPC_metric",
    default = "euclidean",
    type = "character",
    help = "Metric to be used for calculating dissimilarities between observations, available 'euclidean' or 'manhattan' [default : '%default' ]"
  ),
  make_option(
    "--HCPC_method",
    default = "ward",
    type = "character",
    help = "Clustering method between 'ward','average','single', 'complete', 'weighted'  [default :'%default']"
  ),
  make_option(
    "--pdf_out",
    default = "out.pdf",
    type = "character",
    help = "pdf of plots [default : '%default' ]"
  ),
  make_option(
    "--HCPC_consol",
    default = "TRUE",
    type = "logical",
    help = "If TRUE, a k-means consolidation is performed [default :'%default']"
  ),
  make_option(
    "--HCPC_itermax",
    default = "10",
    type = "integer",
    help = "The maximum number of iterations for the consolidation [default :'%default']"
  ),
  make_option(
    "--HCPC_min",
    default = "3",
    type = "integer",
    help = "The least possible number of clusters suggested [default :'%default']"
  ),
  make_option(
    "--HCPC_max",
    default = -1,
    type = "integer",
    help = "The higher possible number of clusters suggested [default :'%default']"
  ),
  make_option(
    "--HCPC_clusterCA",
    default = "rows",
    type = "character",
    help = "A string equals to 'rows' or 'columns' for the clustering of Correspondence Analysis results [default :'%default']"
  ),
  make_option(
    "--HCPC_kk",
    default = Inf,
    type = "numeric",
    help = "The maximum number of iterations for the consolidation [default :'%default']"
  ),
  make_option(
    "--HCPC_clust",
    default = "",
    type = "character",
    help = "Output result of HCPC clustering : two column table (cell identifiers and clusters) [default :'%default']"
  ),
  make_option(
    "--HCPC_mutual_info",
    default = "",
    type = "character",
    help = "Output file of external validation of HCPC clustering with factor levels [default :'%default']"
  ),
  make_option(
    "--HCPC_cluster_description",
    default = "",
    type = "character",
    help = "Output file with variables most contributing to clustering [default :'%default']"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                  args = commandArgs(trailingOnly = TRUE))

if (opt$HCPC_max == -1) {
  opt$HCPC_max <- NULL
}
if (opt$HCPC_kk == -1) {
  opt$HCPC_kk <- Inf
}
#### End of argument implementation ####

#### Input data ####

#We treat data once, at the beginning of the script ####
data <- read.delim(
  opt$data,
  check.names = FALSE,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
# we transpose immediately, because this is the right data structure
data <- as.data.frame(t(data))
#### End of input data ####

######### make PCA with FactoMineR #################
if (opt$visu_choice == "PCA") {
  pdf(opt$pdf_out)
  if (opt$labels) {
    labels <- "ind"
  } else {
    labels <- "none"
  }
  if (opt$factor != "") {
    contrasting_factor <- read.delim(opt$factor, header = TRUE)
    rownames(contrasting_factor) <- contrasting_factor[, 1]
    # we pick only the relevant values of the contrasting factor
    contrasting_factor <- contrasting_factor[rownames(data), ]
    sup <- colnames(contrasting_factor)[2]
    data <- cbind(data, contrasting_factor[, 2])
    colnames(data)[length(data)] <- sup
    if (is.numeric(contrasting_factor[, 2])) {
      res_pca <- PCA(X = data, quanti.sup = sup, graph = FALSE)
      pca_plot <- plot(res_pca, habillage = sup, label = labels,
           title = "PCA graph of cells", cex = opt$item_size,
           axes = c(opt$x_axis, opt$y_axis))
    } else {
      res_pca <- PCA(X = data, quali.sup = sup, graph = FALSE)
      pca_plot <- plot(res_pca, habillage = sup, label = labels,
                       title = "PCA graph of cells", cex = opt$item_size,
                       axes = c(opt$x_axis, opt$y_axis))
    }
  } else {
    res_pca <- PCA(X = data, graph = FALSE)
    pca_plot <- plot(res_pca, label = labels,
                     title = "PCA graph of cells", cex = opt$item_size,
                     axes = c(opt$x_axis, opt$y_axis), col.ind = "deepskyblue4")
  }
  print(pca_plot)
  dev.off()
}
######### END PCA with FactoMineR #################

########### make HCPC with FactoMineR ##########
if (opt$visu_choice == "HCPC") {

# HCPC starts with a PCA
pca <- PCA(
    t(data),
    ncp = opt$HCPC_npc,
    graph = FALSE,
)

pca_ind_coord <- as.data.frame(pca$ind$coord) # coordinates of observations in PCA

# Hierarchical Clustering On Principal Components Followed By Kmean Clustering
res_hcpc <- HCPC(pca,
                 nb.clust = opt$HCPC_ncluster, metric = opt$HCPC_metric, method = opt$HCPC_method,
                 graph = FALSE, consol = opt$HCPC_consol, iter.max = opt$HCPC_itermax, min = opt$HCPC_min, max = opt$HCPC_max,
                 cluster.CA = opt$HCPC_clusterCA, kk = opt$HCPC_kk)
# HCPC plots
dims <- head(as.data.frame(res_hcpc$call$t$res$eig), 2) # dims variances in column 2
pdf(opt$pdf_out)
plot(res_hcpc, choice = "tree")
plot(res_hcpc, choice = "bar")
plot(res_hcpc, choice = "3D.map")
if (opt$labels == FALSE) {
plot(res_hcpc, choice = "map", label = "none")
} else {
plot(res_hcpc, choice = "map")
}

# user contrasts on the pca
if (opt$factor != "") {
  plot(pca, label = "none", col.ind = factor_cols)
  if (is.factor(contrasting_factor$factor)) {
    legend(x = "topright",
       legend = as.character(factor_colors$factor),
       col = factor_colors$color, pch = 16, bty = "n", xjust = 1, cex = 0.7)

    ## Normalized Mutual Information
    sink(opt$HCPC_mutual_info)
    res <- external_validation(
       true_labels = as.numeric(contrasting_factor$factor),
       clusters = as.numeric(res_hcpc$data.clust$clust),
       summary_stats = TRUE
    )
    sink()

  } else {
    legend_col(col = rev(brewer.pal(n = 11, name = "RdYlGn")), lev = cut(contrasting_factor$factor, 11, label = FALSE))
  }
}

dev.off()

if (opt$HCPC_clust != "") {
res_clustering <- data.frame(Cell = rownames(res_hcpc$data.clust),
                             Cluster = res_hcpc$data.clust$clust)

 }

# Description of cluster by most contributing variables / gene expressions

# first transform list of vectors in a list of dataframes
extract_description <- lapply(res_hcpc$desc.var$quanti, as.data.frame)
# second, transfer rownames (genes) to column in the dataframe, before rbinding
extract_description_w_genes <- Map(cbind,
                                   extract_description,
                                   genes = lapply(extract_description, rownames)
                                   )
# Then use data.table to collapse all generated dataframe, with the cluster id in first column
# using the {data.table} rbindlist function
cluster_description <- rbindlist(extract_description_w_genes, idcol = "cluster_id")
cluster_description <- cluster_description[, c(8, 1, 2, 3, 4, 5, 6, 7)] # reorganize columns


# Finally, output cluster description data frame
write.table(
    cluster_description,
    file = opt$HCPC_cluster_description,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
    )

}

## Return cluster table to user

if (opt$HCPC_clust != "") {
  write.table(
    res_clustering,
    file = opt$HCPC_clust,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
    )
}

################  t-SNE ####################
if (opt$visu_choice == "tSNE") {
  # filter and transpose df for tsne and pca
  tdf <- t(data)
  # make tsne and plot results
  set.seed(opt$Rtsne_seed) ## Sets seed for reproducibility

  tsne_out <- Rtsne(tdf,
                    dims = opt$Rtsne_dims,
                    initial_dims = opt$Rtsne_initial_dims,
                    perplexity = opt$Rtsne_perplexity,
                    theta = opt$Rtsne_theta,
                    max_iter = opt$Rtsne_max_iter,
                    pca = opt$Rtsne_pca,
                    pca_center = opt$Rtsne_pca_center,
                    pca_scale = opt$Rtsne_pca_scale,
                    normalize = opt$Rtsne_normalize,
                    exaggeration_factor = opt$Rtsne_exaggeration_factor)

  embedding <- as.data.frame(tsne_out$Y[, 1:2])
  embedding$Class <- as.factor(rownames(tdf))
  gg_legend <- theme(legend.position = "right")
  if (opt$factor == "") {
    ggplot(embedding, aes(x = V1, y = V2)) +
      geom_point(size = 1, color = "deepskyblue4") +
      gg_legend +
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      ggtitle("t-SNE") +
      if (opt$labels) {
        geom_text(aes(label = Class), hjust = -0.2, vjust = -0.5, size = 1.5, color = "deepskyblue4")
      }
    } else {
    if (is.numeric(contrasting_factor$factor)) {
      embedding$factor <- contrasting_factor$factor
    } else {
      embedding$factor <- as.factor(contrasting_factor$factor)
    }

    ggplot(embedding, aes(x = V1, y = V2, color = factor)) +
      geom_point(size = 1) +
      gg_legend +
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      ggtitle("t-SNE") +
      if (opt$labels) {
        geom_text(aes(label = Class, colour = factor), hjust = -0.2, vjust = -0.5, size = 1.5)
      }
    }
  ggsave(file = opt$pdf_out, device = "pdf")

}

