##################################################################################################
# Running PATHIFIER (Drier et al., 2013)
# Based on the work of Author: Miguel Angel Garcia-Campos - Github: https://github.com/AngelCampos
##################################################################################################


options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
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
    help = "File separator [default : '%default' ]"
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
    help = "Log file name [default : '%default' ]"
  ),
  make_option(
    "--max_stability",
    type = "logical",
    default = T,
    help = "If true, throw away components leading to low stability of sampling noise [default : '%default' ]"
  ),
  make_option(
    "--attempts",
    type = "integer",
    default = 10,
    help = "Number of runs to determine stability. [default : '%default' ]"
  ),
  make_option(
    "--min_std",
    type = "character",
    default = "0.4",
    help = "Minimum of standard deviation to filter out low variable genes. 
    Use --min.std data, to use the minimum std of your data [default : '%default' ]"
  ),
  make_option(
    "--min_exp",
    type = "character",
    default = "4",
    help = "Minimum of gene expression to filter out low expressed genes. 
    Use --min.exp data, to use the minimum expression of your data [default : '%default' ]"
  ),
  make_option(
    "--pds",
    type = "character",
    default = "PDS.tsv",
    help = "Output PDS (Pathway deregulation score) of Pathifier in tabular file [default : '%default' ]"
  ),
  make_option(
    "--heatmap_cluster_cells",
    type = "logical",
    default = TRUE,
    help = "Cluster columns (cells) in the heatmap [default : '%default' ]"
  ),
  make_option(
    "--heatmap_cluster_pathways",
    type = "logical",
    default = TRUE,
    help = "Cluster rows (pathways) in the heatmap [default : '%default' ]"
  ),
  make_option(
    "--heatmap_show_cell_labels",
    type = "logical",
    default = FALSE,
    help = "Print column names (cells) on the heatmap [default : '%default' ]"
  ),
  make_option(
    "--heatmap_show_pathway_labels",
    type = "logical",
    default = FALSE,
    help = "Print row names (pathways) on the heatmap [default : '%default' ]"
  ),
  make_option(
    "--plot",
    type = "character",
    default = "./plot.pdf",
    help = "Pathifier visualization [default : '%default' ]"
  )
)
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)
if (args$sep == "tab") {args$sep = "\t"}


# set seed for reproducibility
set.seed(123)

# Load expression data for PATHIFIER
exp.matrix <- as.matrix(read.delim(file = args$exp, as.is = T, row.names = 1, sep = args$sep, check.names = F))

# Load Genesets annotation 
gene_sets_file <- file(args$genes,open="r")
gene_sets <- readLines(gene_sets_file)
close(gene_sets_file)

#  Generate a list that contains genes in genesets
gs <- strsplit(gene_sets, "\t")
names(gs) <- lapply(gs, function(x) x[1])
gs <- lapply(gs, function(x) x[-c(1:2)])

# Generate a list that contains the names of the genesets used
# pathwaynames <- as.list(gene_sets[,1])
pathwaynames <- names(gs)

# Prepare data and parameters ##################################################
# Extract information from binary phenotypes. 1 = Normal, 0 = Tumor
if(args$is_normal == T){
  normals <- read.delim(file = args$normals, h = F)
  normals <- as.logical(normals[,2])
  N.exp.matrix <- exp.matrix[,normals]
} else N.exp.matrix <- exp.matrix

# Calculate MIN_STD
rsd <- apply(N.exp.matrix, 1, sd)
min_std <- quantile(rsd, 0.25) 

# Calculate MIN_EXP 
min_exp <- quantile(as.vector(as.matrix(exp.matrix)), 0.1) # Percentile 10 of data

# Filter low value genes. At least 10% of samples with values over min_exp
# Set expression levels < MIN_EXP to MIN_EXP
over <- apply(exp.matrix, 1, function(x) x > min_exp)
G.over <- apply(over, 2, mean)
G.over <- names(G.over)[G.over > 0.1]
exp.matrix.filtered <- exp.matrix[G.over,]
exp.matrix.filtered[exp.matrix.filtered < min_exp] <- min_exp

# Set maximum 5000 genes with more variance
V <- names(sort(apply(exp.matrix.filtered, 1, var), decreasing = T))[1:5000]
V <- V[!is.na(V)]
exp.matrix.filtered <- exp.matrix.filtered[V,]
allgenes <- as.vector(rownames(exp.matrix.filtered))


if(args$min_std == "data"){
  args$min_std <- min_std
} else args$min_std <- as.numeric(args$min_std)

if(args$min_exp == "data"){
  args$min_exp <- min_exp
} else args$min_exp <- as.numeric(args$min_exp)


# Open pdf
pdf(args$plot)

# Construct continuous color scale
myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

col_score_fun = colorRamp2(c(0, 0.5, 1), c("blue", "yellow", "red"))

#Continous color scale legend
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

# Run Pathifier
if(args$is_normal == T){
  PDS <- quantify_pathways_deregulation(exp.matrix.filtered, 
                                        allgenes,
                                        gs,
                                        pathwaynames,
                                        normals, 
                                        maximize_stability = args$max_stability,
                                        attempts = args$attempts,
                                        logfile = args$logfile,
                                        min_std = args$min_std,
                                        min_exp = args$min_exp)

  for(i in pathwaynames){
    DF <- data.frame(PDS$curves[[i]][,1:3], normal = normals, PDS = as.numeric(PDS$scores[[i]]), curve_order = as.vector(PDS$curves_order[[i]]))
    ordered <- DF[DF$curve_order,]

    cols_score <- myColorRamp(c("blue", "red"), ordered$PDS)
    sc3_score <- scatterplot3d(ordered[,1:3], main = paste("Principal curve of", i), box = F, pch = 19, type = "l")
    sc3_score$points3d(ordered[,1:3], box = F, pch = 19, col = cols_score)
    #legend.col(cols_score, levels = ordered$PDS)

    cols_status <- ifelse(ordered$normal, "blue", "red")
    sc3 <- scatterplot3d(ordered[,1:3], main = paste("Principal curve of", i), box = F, pch = 19, type = "l")
    sc3$points3d(ordered[,1:3], box = F, pch = ifelse(ordered$normal, 19, 8), col = cols_status)
  }

} else{
  PDS <- quantify_pathways_deregulation(exp.matrix.filtered, 
                                        allgenes,
                                        gs,
                                        pathwaynames,
                                        maximize_stability = args$max_stability,
                                        attempts = args$attempts,
                                        logfile = args$logfile,
                                        min_std = args$min_std,
                                        min_exp = args$min_exp)
  
  for(i in pathwaynames){
    DF <- data.frame(PDS$curves[[i]][,1:3], PDS = as.numeric(PDS$scores[[i]]), curve_order = as.vector(PDS$curves_order[[i]]))
    ordered <- DF[DF$curve_order,]

    cols <- myColorRamp(c("blue", "yellow", "red"), ordered$PDS) 

    sc3 <- scatterplot3d(ordered[,1:3], main = paste("Principal curve of", i), box = F, pch = 19, type = "l")
    sc3$points3d(ordered[,1:3], box = F, pch = 19, col = col_score_fun(ordered$PDS))
    #legend.col(col_score_fun(ordered$PDS), levels = ordered$PDS)
    sc3.coords <- sc3$xyz.convert(ordered[,1:3])
    text(sc3.coords$x, sc3.coords$y, labels = as.character(round(ordered$PDS, 3)), pos = 2, offset = 0.5, box = F)

    sc3 <- scatterplot3d(ordered[,1:3], main = paste("Principal curve of", i), box = F, pch = 19, type = "l")
    sc3$points3d(ordered[,1:3], box = F, pch = 19, col = cols)
    #legend.col(cols, levels = ordered$PDS)
    sc3.coords <- sc3$xyz.convert(ordered[,1:3])
    text(sc3.coords$x, sc3.coords$y, labels = as.character(round(ordered$PDS, 3)), pos = 2, offset = 0.5, box = F)

  }
}


PDS_scores <- mapply(FUN = cbind, PDS$scores)
dimnames(PDS_scores) <- list(colnames(exp.matrix.filtered), names(PDS$scores))



## plot heatmap 
if(ncol(PDS_scores) > 1){
    pheatmap(t(PDS_scores),
             main = "Heatmap of Pathway Deregulation Score",   # heat map title
             cluster_rows = args$heatmap_cluster_pathways,     # apply clustering method
             cluster_cols = args$heatmap_cluster_cells,        # apply clustering method
                   
             #Additional Options
             ## Color labeled columns
             # annotation_col = cell_metadata,
             # annotation_colors = cell_metadata_color,
             show_rownames = args$heatmap_show_pathway_labels,
             show_colnames = args$heatmap_show_cell_labels,
             border_color = NA,
             legend = TRUE
    )
}
dev.off()


## write table
write.table(PDS_scores,
            args$pds,
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t"
)
