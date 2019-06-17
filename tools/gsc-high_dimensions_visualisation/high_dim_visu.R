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
    help = "add labels to t-SNE plot [default : '%default' ]"
  ),
  make_option(
    "--visu_choice",
    default = 'PCA',
    type = 'character',
    help = "visualisation method ('PCA', 'tSNE', 'HCPC') [default : '%default' ]"
  ),
  make_option(
    "--plot_coordinates",
    default = FALSE,
    type = 'logical',
    help = "Table with plot coordinates in the output [default : '%default' ]"
  ),
   make_option(
    "--table_coordinates",
    default = 'Coord.tab',
    type = 'character',
    help = "Table with plot coordinates [default : '%default' ]"
  ),
  make_option(
    "--seed",
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
    "--PCA_ncp",
    default = 5,
    type = 'integer',
    help = "number of dimensions kept in the results [default : '%default' ]"
  ),
  make_option(
    "--HCPC_ncluster",
    default = -1,
    type = 'numeric',
    help = "nb.clust, number of clusters to consider in the hierarchical clustering. [default : -1 let HCPC to optimize the number]"
  ),
   make_option(
    "--HCPC_ncp",
    default = 5,
    type = 'integer',
    help = "npc, number of dimensions which are kept for HCPC analysis [default : '%default' ]"
  ),
  make_option(
    "--HCPC_metric",
    default = 'euclidian',
    type = 'character',
    help = "Metric to be used for calculating dissimilarities between observations, available 'euclidian' or 'manhattan' [default : '%default' ]"
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
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}
if (opt$sep == "comma") {opt$sep = ","}

data = read.table(
  opt$data,
  check.names = FALSE,
  header = opt$colnames,
  row.names = 1,
  sep = opt$sep
)

# t-SNE
if (opt$visu_choice == 'tSNE') {
  # filter and transpose df for tsne and pca
  data = data[rowSums(data) != 0,] # remove lines without information (with only 0s)
  tdf = t(data)
  # make tsne and plot results
  set.seed(opt$seed) ## Sets seed for reproducibility

  tsne_out <- Rtsne(tdf, dims = opt$Rtsne_dims, initial_dims = opt$Rtsne_initial_dims, 
     perplexity = opt$Rtsne_perplexity , theta = opt$Rtsne_theta, pca = opt$Rtsne_pca, 
     pca_center = opt$Rtsne_pca_center,  pca_scale = opt$Rtsne_pca_scale,
     normalize = opt$Rtsne_normalize, exaggeration_factor=opt$Rtsne_exaggeration_factor)

  embedding <- as.data.frame(tsne_out$Y[,1:2])
  embedding$Class <- as.factor(rownames(tdf))
  gg_legend = theme(legend.position="none")
  ggplot(embedding, aes(x=V1, y=V2)) +
    geom_point(size=1, color='red') +
    gg_legend +
    xlab("") +
    ylab("") +
    ggtitle('t-SNE') +
    if (opt$labels == TRUE) {
      geom_text(aes(label=Class),hjust=-0.2, vjust=-0.5, size=1.5, color='darkblue')
    }
  ggsave(file=opt$pdf_out, device="pdf")
  
  #save coordinates table
  if(opt$plot_coordinates){
  coord_table <- cbind(rownames(tdf),as.data.frame(tsne_out$Y))
  colnames(coord_table)=c("Cells",paste0("DIM",(1:opt$Rtsne_dims)))
  }
}


# make PCA with FactoMineR
if (opt$visu_choice == 'PCA') {
  pca <- PCA(t(data), ncp=opt$PCA_ncp, graph=FALSE)
  pdf(opt$pdf_out)
  if (opt$labels == FALSE) {
    plot(pca, label="none")
    } else {
    plot(pca, cex=0.2)
  }
dev.off()

  #save coordinates table
  if(opt$plot_coordinates){
  coord_table <- cbind(rownames(pca$ind$coord),as.data.frame(pca$ind$coord))
  colnames(coord_table)=c("Cells",paste0("DIM",(1:opt$PCA_ncp)))
  }

}

#############################
# make HCPC with FactoMineR #
if (opt$visu_choice == 'HCPC') {
pdf(opt$pdf_out)

## HCPC starts with a PCA
pca <- PCA(
    t(data),
    ncp = opt$HCPC_ncp,
    graph = FALSE,
    scale.unit = FALSE
)

PCA_IndCoord = as.data.frame(pca$ind$coord) # coordinates of observations in PCA
#Vartab = get_eig(pca)  #variable masquee car on ne s en sert pas dans la suite du code




## Hierarchical Clustering On Principal Components Followed By Kmean Clustering
res.hcpc <- HCPC(pca,
                 nb.clust=opt$HCPC_ncluster, metric=opt$HCPC_metric, method=opt$HCPC_method,
                 graph=F)
# A string. "tree" plots the tree. "bar" plots bars of inertia gains. "map" plots a factor map,
# individuals colored by cluster. "3D.map" plots the same factor map, individuals colored by cluster,
# the tree above.
plot(res.hcpc, choice="tree")
plot(res.hcpc, choice="bar")
plot(res.hcpc, choice="3D.map")
if (opt$labels == FALSE) {
plot(res.hcpc, choice="map", label="none")
} else {
plot(res.hcpc, choice="map")
}


QuanVarDescCluster <- as.list(res.hcpc$desc.var$quanti)
AxesDescCluster = as.list(res.hcpc$desc.axes$quanti)


## Clusters to which individual observations belong
Clust <- data.frame(Cluster = res.hcpc$data.clust[, (nrow(data) + 1)],
                    Observation = rownames(res.hcpc$data.clust))
metadata <- data.frame(Observation=colnames(data), row.names=colnames(data))
metadata = merge(y = metadata,
                 x = Clust,
                 by = "Observation")
ObsNumberPerCluster = as.data.frame(table(metadata$Cluster))
colnames(ObsNumberPerCluster) = c("Cluster", "ObsNumber")

## Silhouette Plot
hc.cut = hcut(PCA_IndCoord,
              k = nlevels(metadata$Cluster),
              hc_method = "ward.D2")
              
Sil = fviz_silhouette(hc.cut)
sil1 = as.data.frame(Sil$data)

## Normalized Mutual Information # to be implemented later
# sink(opt$mutual_info)
# res = external_validation(
#   as.numeric(factor(metadata[, Patient])),
#   as.numeric(metadata$Cluster),
#   method = "adjusted_rand_index",
#   summary_stats = TRUE
# )
# sink()

#  plot(pca, label="none", habillage="ind", col.hab=metadata$Cluster_2d_color )
#  plot(pca, label="none", habillage="ind", col.hab=cols )
#  PatientSampleColor = Color1[as.factor(metadata[, Patient])]
#  plot(pca, label="none", habillage="ind", col.hab=PatientSampleColor )
dev.off()

 if(opt$plot_coordinates){
  coord_table <- cbind(Cell=rownames(res.hcpc$call$X),as.data.frame(res.hcpc$call$X))
  colnames(coord_table)=c("Cells",paste0("DIM",(1:opt$HCPC_ncp)),"Cluster")
  }




}

## Return coordinates file to user

if(opt$plot_coordinates){
  write.table(
    coord_table,
    file = opt$table_coordinates,
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
    )
}


  











