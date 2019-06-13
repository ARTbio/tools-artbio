if (length(commandArgs(TRUE)) == 0) {
  system("Rscript cpm_tpm_rpk.R -h", intern = F)
  q("no")
}


# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)
library(ggplot2)
library(reshape2)
library(Rtsne)
library(ggfortify)



#Arguments
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
    "--tsne_labels",
    default = FALSE,
    type = 'logical',
    help = "add labels to t-SNE plot [default : '%default' ]"
  ),
  make_option(
    "--seed",
    default = 42,
    type = 'integer',
    help = "Seed value for reproducibility [default : '%default' ]"
  ),
  make_option(
    "--perp",
    default = 5.0,
    type = 'numeric',
    help = "perplexity [default : '%default' ]"
  ),
  make_option(
    "--theta",
    default = 1.0,
    type = 'numeric',
    help = "theta [default : '%default' ]"
  ),
  make_option(
    "--tsne_out",
    default = "tsne.pdf",
    type = 'character',
    help = "T-SNE pdf [default : '%default' ]"
  ),
  make_option(
    "--pca_out",
    default = "pca.pdf",
    type = 'character',
    help = "PCA pdf [default : '%default' ]"
  )

)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$sep == "tab") {opt$sep = "\t"}

data = read.table(
  opt$data,
  check.names = FALSE,
  header = opt$colnames,
  row.names = 1,
  sep = opt$sep
)

colnames(res) = colnames(data)


if (opt$log == TRUE) {
  res = log2(res + 1)
}

write.table(
  cbind(Features = rownames(res), res),
  opt$out,
  col.names = opt$colnames,
  row.names = F,
  quote = F,
  sep = "\t"
)

## 
if (opt$visu == TRUE) {
  df = res
  # filter and transpose df for tsne and pca
  df = df[rowSums(df) != 0,] # remove lines without information (with only 0 counts)
  tdf = t(df)
  # make tsne and plot results
  set.seed(opt$seed) ## Sets seed for reproducibility
  tsne_out <- Rtsne(tdf, perplexity=opt$perp, theta=opt$theta) # 
  embedding <- as.data.frame(tsne_out$Y)
  embedding$Class <- as.factor(sub("Class_", "", rownames(tdf)))
  gg_legend = theme(legend.position="none")
  ggplot(embedding, aes(x=V1, y=V2)) +
    geom_point(size=1, color='red') +
    gg_legend +
    xlab("") +
    ylab("") +
    ggtitle('t-SNE') +
    if (opt$tsne_labels == TRUE) {
      geom_text(aes(label=Class),hjust=-0.2, vjust=-0.5, size=2.5, color='darkblue')
    }
  ggsave(file=opt$tsne_out, device="pdf")
  # make PCA and plot result with ggfortify (autoplot)
  tdf.pca <- prcomp(tdf, center = TRUE, scale. = T)
  if (opt$tsne_labels == TRUE) {
      autoplot(tdf.pca, shape=F, label=T, label.size=2.5, label.vjust=1.2,
               label.hjust=1.2,
               colour="darkblue") +
      geom_point(size=1, color='red') +
      xlab(paste("PC1",summary(tdf.pca)$importance[2,1]*100, "%")) +
      ylab(paste("PC2",summary(tdf.pca)$importance[2,2]*100, "%")) +
      ggtitle('PCA')
      ggsave(file=opt$pca_out, device="pdf")   
      } else {
      autoplot(tdf.pca, shape=T, colour="darkblue") +
      geom_point(size=1, color='red') +
      xlab(paste("PC1",summary(tdf.pca)$importance[2,1]*100, "%")) +
      ylab(paste("PC2",summary(tdf.pca)$importance[2,2]*100, "%")) +
      ggtitle('PCA') 
      ggsave(file=opt$pca_out, device="pdf")
  }
}
  











