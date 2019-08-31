# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
requiredPackages = c('optparse', 'Rtsne', 'ggplot2', 'ggfortify')
warnings()
library(optparse)
library(FactoMineR)
library(factoextra)
library(Rtsne)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(ClusterR)

# Arguments
option_list = list(
  make_option(
    "--data",
    default = NA,
    type = 'character',
    help = "Input file that contains expression value to visualise"
  ),
  make_option(
    "--sep",
    default = '\t',
    type = 'character',
    help = "File separator [default : '%default' ]"
  ),
  make_option(
    "--colnames",
    default = TRUE,
    type = 'logical',
    help = "Consider first line as header ? [default : '%default' ]"
  ),
  make_option(
    "--out",
    default = "res.tab",
    type = 'character',
    help = "Output name [default : '%default' ]"
  ),
  make_option(
    "--labels",
    default = FALSE,
    type = 'logical',
    help = "add labels in scatter plots [default : '%default' ]"
  ),
  make_option(
    "--factor",
    default = '',
    type = 'character',
    help = "A two column table that specifies factor levels for contrasting data [default : '%default' ]"
  ),
  make_option(
    "--visu_choice",
    default = 'PCA',
    type = 'character',
    help = "visualisation method ('PCA', 'tSNE', 'HCPC') [default : '%default' ]"
  ),
  make_option(
    "--table_coordinates",
    default = '',
    type = 'character',
    help = "Table with plot coordinates [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_seed",
    default = 42,
    type = 'integer',
    help = "Seed value for reproducibility [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_dims",
    default = 2,
    type = 'integer',
    help = "Output dimensionality [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_initial_dims",
    default = 50,
    type = 'integer',
    help = "The number of dimensions that should be retained in the initial PCA step [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_perplexity",
    default = 5.0,
    type = 'numeric',
    help = "perplexity [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_theta",
    default = 1.0,
    type = 'numeric',
    help = "theta [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_max_iter",
    default = 1000,
    type = 'integer',
    help = "max_iter [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_pca",
    default = TRUE,
    type = 'logical',
    help = "Whether an initial PCA step should be performed [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_pca_center",
    default = TRUE,
    type = 'logical',
    help = "Should data be centered before pca is applied? [default : '%default' ]"
  ),
   make_option(
    "--Rtsne_pca_scale",
    default = FALSE,
    type = 'logical',
    help = "Should data be scaled before pca is applied? [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_normalize",
    default = TRUE,
    type = 'logical',
    help = "Should data be normalized internally prior to distance calculations? [default : '%default' ]"
  ),
  make_option(
    "--Rtsne_exaggeration_factor",
    default = 12.0,
    type = 'numeric',
    help = " Exaggeration factor used to multiply the P matrix in the first part of the optimization [default : '%default' ]"
  ),
   make_option(
    "--PCA_npc",
    default = 5,
    type = 'integer',
    help = "number of dimensions kept in the results [default : '%default' ]"
  ),
   make_option(
    "--PCA_x_axis",
    default = 1,
    type = 'integer',
    help = "PC to plot in the x axis [default : '%default' ]"
  ),
   make_option(
    "--PCA_y_axis",
    default = 2,
    type = 'integer',
    help = "PC to plot in the y axis [default : '%default' ]"
  ),
  make_option(
    "--HCPC_ncluster",
    default = -1,
    type = 'numeric',
    help = "nb.clust, number of clusters to consider in the hierarchical clustering. [default : -1 let HCPC to optimize the number]"
  ),
   make_option(
    "--HCPC_npc",
    default = 5,
    type = 'integer',
    help = "npc, number of dimensions which are kept for HCPC analysis [default : '%default' ]"
  ),
  make_option(
    "--HCPC_metric",
    default = 'euclidean',
    type = 'character',
    help = "Metric to be used for calculating dissimilarities between observations, available 'euclidean' or 'manhattan' [default : '%default' ]"
  ),
  make_option(
    "--HCPC_method",
    default = 'ward',
    type = 'character',
    help = "Clustering method between 'ward','average','single', 'complete', 'weighted'  [default :'%default']"
  ),
  make_option(
    "--pdf_out",
    default = "out.pdf",
    type = 'character',
    help = "pdf of plots [default : '%default' ]"
  ),
  make_option(
    "--HCPC_consol",
    default = 'TRUE',
    type = 'logical',
    help = "If TRUE, a k-means consolidation is performed [default :'%default']"
  ),
  make_option(
    "--HCPC_itermax",
    default = '10',
    type = 'integer',
    help = "The maximum number of iterations for the consolidation [default :'%default']"
  ),
  make_option(
    "--HCPC_min",
    default = '3',
    type = 'integer',
    help = "The least possible number of clusters suggested [default :'%default']"
  ),
   make_option(
    "--HCPC_max",
    default = -1,
    type = 'integer',
    help = "The higher possible number of clusters suggested [default :'%default']"
  ),
   make_option(
    "--HCPC_clusterCA",
    default = 'rows',
    type = 'character',
    help = "A string equals to 'rows' or 'columns' for the clustering of Correspondence Analysis results [default :'%default']"
  ),
  make_option(
    "--HCPC_kk",
    default = Inf,
    type = 'numeric',
    help = "The maximum number of iterations for the consolidation [default :'%default']"
  ),
  make_option(
    "--HCPC_clust",
    default = "",
    type = 'character',
    help = "Output result of HCPC clustering : two column table (cell identifiers and clusters) [default :'%default']"
  ),
  make_option(
    "--mutual_info",
    default = "",
    type = "character",
    help = "Output file of external validation of HCPC clustering with factor levels [default :'%default']"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep <- "\t"}
if (opt$sep == "comma") {opt$sep <- ","}
if(opt$HCPC_max == -1) {opt$HCPC_max <- NULL}
if(opt$HCPC_kk == -1) {opt$HCPC_kk <- Inf}

##Add legend to plot()
legend.col <- function(col, lev){

opar <- par

n <- length(col)

bx <- par("usr")

box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
box.cy <- c(bx[3], bx[3])
box.sy <- (bx[4] - bx[3]) / n

xx <- rep(box.cx, each = 2)

par(xpd = TRUE)
for(i in 1:n){

yy <- c(box.cy[1] + (box.sy * (i - 1)),
box.cy[1] + (box.sy * (i)),
box.cy[1] + (box.sy * (i)),
box.cy[1] + (box.sy * (i - 1)))
polygon(xx, yy, col = col[i], border = col[i])

}
par(new = TRUE)
plot(0, 0, type = "n",
ylim = c(min(lev), max(lev)),
yaxt = "n", ylab = "",
xaxt = "n", xlab = "",
frame.plot = FALSE)
axis(side = 4, las = 2, tick = FALSE, line = .25)
par <- opar
}


data = read.table(
  opt$data,
  check.names = FALSE,
  header = opt$colnames,
  row.names = 1,
  sep = opt$sep
)

# Contrasting factor and its colors
if (opt$factor != '') {
  contrasting_factor <- read.delim(
    opt$factor,
    header = TRUE
  )
  rownames(contrasting_factor) <- contrasting_factor[,1]
  contrasting_factor <- contrasting_factor[colnames(data),]
  colnames(contrasting_factor) <- c("id","factor")
  if(is.numeric(contrasting_factor$factor)){
    factor_cols <- rev(brewer.pal(n = 11, name = "RdYlGn"))[contrasting_factor$factor]
  } else {
    contrasting_factor$factor <- as.factor(contrasting_factor$factor)
    if(nlevels(contrasting_factor$factor) == 2){
      colors_labels <- c("#E41A1C", "#377EB8")
    } else {
      colors_labels <- brewer.pal(nlevels(contrasting_factor$factor), name = 'Paired')
    }
    factorColors <-
      with(contrasting_factor,
           data.frame(factor = levels(contrasting_factor$factor),
                      color = I(colors_labels)
           )
      )
    factor_cols <- factorColors$color[match(contrasting_factor$factor,
                                          factorColors$factor)]
  }
} else {
  factor_cols <- rep("deepskyblue4", length(rownames(data)))
}

################  t-SNE ####################
if (opt$visu_choice == 'tSNE') {
  # filter and transpose df for tsne and pca
  tdf = t(data)
  # make tsne and plot results
  set.seed(opt$Rtsne_seed) ## Sets seed for reproducibility

  tsne_out <- Rtsne(tdf,
                    dims = opt$Rtsne_dims,
                    initial_dims = opt$Rtsne_initial_dims, 
                    perplexity = opt$Rtsne_perplexity ,
                    theta = opt$Rtsne_theta,
                    max_iter = opt$Rtsne_max_iter,
                    pca = opt$Rtsne_pca, 
                    pca_center = opt$Rtsne_pca_center,
                    pca_scale = opt$Rtsne_pca_scale,
                    normalize = opt$Rtsne_normalize,
                    exaggeration_factor=opt$Rtsne_exaggeration_factor)

  embedding <- as.data.frame(tsne_out$Y[,1:2])
  embedding$Class <- as.factor(rownames(tdf))
  gg_legend = theme(legend.position="right")
  if (opt$factor == '') {
    ggplot(embedding, aes(x=V1, y=V2)) +
      geom_point(size=1, color='deepskyblue4') +
      gg_legend +
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      ggtitle('t-SNE') +
      if (opt$labels) {
        geom_text(aes(label=Class),hjust=-0.2, vjust=-0.5, size=1.5, color='deepskyblue4')
      }
    } else {
    if(is.numeric(contrasting_factor$factor)){
      embedding$factor <- contrasting_factor$factor
    } else {
      embedding$factor <- as.factor(contrasting_factor$factor)
    }

    ggplot(embedding, aes(x=V1, y=V2, color=factor)) +
      geom_point(size=1) + #, color=factor_cols
      gg_legend +
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      ggtitle('t-SNE') +
      if (opt$labels) {
        geom_text(aes(label=Class, colour=factor),hjust=-0.2, vjust=-0.5, size=1.5)
      }
    }    
  ggsave(file=opt$pdf_out, device="pdf")
  
  #save coordinates table
  if(opt$table_coordinates != ''){
  coord_table <- cbind(rownames(tdf), round(as.data.frame(tsne_out$Y), 6))
  colnames(coord_table)=c("Cells",paste0("DIM",(1:opt$Rtsne_dims)))
  }
}


######### make PCA with FactoMineR #################
if (opt$visu_choice == 'PCA') {
  pca <- PCA(t(data), ncp=opt$PCA_npc, graph=FALSE)
  pdf(opt$pdf_out)
  if (opt$labels == FALSE) {
    plot(pca, axes = c(opt$PCA_x_axis,opt$PCA_y_axis), label="none" , col.ind = factor_cols)
    } else {
    plot(pca, axes = c(opt$PCA_x_axis,opt$PCA_y_axis), cex=0.2 , col.ind = factor_cols)
  }
if (opt$factor != '') {
  if(is.factor(contrasting_factor$factor)) {
    legend(x = 'topright', 
       legend = as.character(factorColors$factor),
       col = factorColors$color, pch = 16, bty = 'n', xjust = 1, cex=0.7)
  } else {
    legend.col(col = rev(brewer.pal(n = 11, name = "RdYlGn")), lev = cut(contrasting_factor$factor, 11, label = FALSE))
  }
}
dev.off()

  #save coordinates table
  if(opt$table_coordinates != ''){
  coord_table <- cbind(rownames(pca$ind$coord), round(as.data.frame(pca$ind$coord), 6))
  colnames(coord_table)=c("Cells",paste0("DIM",(1:opt$PCA_npc)))
  }

}

########### make HCPC with FactoMineR ##########
if (opt$visu_choice == 'HCPC') {

# HCPC starts with a PCA
pca <- PCA(
    t(data),
    ncp = opt$HCPC_npc,
    graph = FALSE,
)

PCA_IndCoord = as.data.frame(pca$ind$coord) # coordinates of observations in PCA

# Hierarchical Clustering On Principal Components Followed By Kmean Clustering
res.hcpc <- HCPC(pca,
                 nb.clust=opt$HCPC_ncluster, metric=opt$HCPC_metric, method=opt$HCPC_method,
                 graph=F,consol=opt$HCPC_consol,iter.max=opt$HCPC_itermax,min=opt$HCPC_min,max=opt$HCPC_max,
                 cluster.CA=opt$HCPC_clusterCA,kk=opt$HCPC_kk)
# HCPC plots
dims <- head(as.data.frame(res.hcpc$call$t$res$eig),2) # dims variances in column 2
pdf(opt$pdf_out)
plot(res.hcpc, choice="tree")
plot(res.hcpc, choice="bar")
plot(res.hcpc, choice="3D.map")
if (opt$labels == FALSE) {
plot(res.hcpc, choice="map", label="none")
} else {
plot(res.hcpc, choice="map")
}

# user contrasts on the pca
if (opt$factor != '') {
  plot(pca, label="none", col.ind = factor_cols)
  if(is.factor(contrasting_factor$factor)) {
    legend(x = 'topright', 
       legend = as.character(factorColors$factor),
       col = factorColors$color, pch = 16, bty = 'n', xjust = 1, cex=0.7)
    
    ## Normalized Mutual Information
    sink(opt$mutual_info)
    res <- external_validation(
       true_labels = as.numeric(contrasting_factor$factor),
       clusters = as.numeric(res.hcpc$data.clust$clust),
       summary_stats = TRUE
    )
    sink()

  } else {
    legend.col(col = rev(brewer.pal(n = 11, name = "RdYlGn")), lev = cut(contrasting_factor$factor, 11, label = FALSE))
  }
}
## Clusters to which individual observations belong # used ?
# Clust <- data.frame(Cluster = res.hcpc$data.clust[, (nrow(data) + 1)],
#                     Observation = rownames(res.hcpc$data.clust))
# metadata <- data.frame(Observation=colnames(data), row.names=colnames(data))
# metadata = merge(y = metadata,
#                  x = Clust,
#                  by = "Observation")

# unclear utility
# ObsNumberPerCluster = as.data.frame(table(metadata$Cluster))
# colnames(ObsNumberPerCluster) = c("Cluster", "ObsNumber")
# 
# ## Silhouette Plot # not used
# hc.cut = hcut(PCA_IndCoord,
#               k = nlevels(metadata$Cluster),
#               hc_method = "ward.D2")
#               
# Sil = fviz_silhouette(hc.cut)
# sil1 = as.data.frame(Sil$data)

dev.off()

 if(opt$table_coordinates != ''){
  coord_table <- cbind(Cell=rownames(res.hcpc$call$X),
                       round(as.data.frame(res.hcpc$call$X[, -length(res.hcpc$call$X)]), 6),
                       as.data.frame(res.hcpc$call$X[, length(res.hcpc$call$X)])
                       )
  colnames(coord_table)=c("Cells",paste0("DIM",(1:opt$HCPC_npc)),"Cluster")
  }

 if(opt$HCPC_clust != ""){
 res_clustering <- data.frame(Cell = rownames(res.hcpc$data.clust),
                              Cluster = res.hcpc$data.clust$clust)
 
 }

}

## Return coordinates file to user

if(opt$table_coordinates != ''){
  write.table(
    coord_table,
    file = opt$table_coordinates,
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
    )
}


if(opt$HCPC_clust != ""){
  write.table(
    res_clustering,
    file = opt$HCPC_clust,
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
    )
}  











