## Setup R error handling to go to stderr
options(
  show.error.messages = F,
  error = function () {
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
  }
)

option_list <- list(
  optparse::make_option(c("-f", "--file"),
                        type = "character",
                        help = "Path to the input file"),
  optparse::make_option(c("-c", "--colnames"),
                        type = "logical",
                        help = "Consider first line as header"),
  optparse::make_option(c("-s", "--sep"),
                        type = "character",
                        help = "Input file separator"),
  optparse::make_option(c("-m", "--is_mito"),
                        type = "logical",
                        help = "There are mitochondrial genes in the dataset"),
  optparse::make_option("--header_mito",
                        type = "character",
                        help = "Prefixe of mitochondrial genes in the dataset.
                        e.g. 'MT-'"),
  optparse::make_option(c("-p", "--package"),
                        type = "character",
                        help = "Name of the package used for the filtering :
                        scater or Seurat"),
  optparse::make_option(c("-g", "--gene_filter"),
                        type = "character",
                        help = "Type of filtering :
                        low.abundances or min.cells"),
  optparse::make_option("--min_cells",
                        type = "double",
                        help = "Threshold for the gene filtering. e.g. '3'"),
  optparse::make_option("--min_genes",
                        type = "double",
                        help = "Threshold for the cell filtering. e.g. '5'"),
  optparse::make_option("--output_matrix",
                        type = "character",
                        help = "Path to the filtered matrix"),
  optparse::make_option("--output_pdf",
                        type = "character",
                        help = "Path to the pdf file with plots")
)

parser <- optparse::OptionParser(usage = "%prog [options] file",
                                 option_list = option_list)
args <- optparse::parse_args(parser)



##Import dataset
data.counts <- read.table(
  args$file,
  header = args$colnames,
  stringsAsFactors = F,
  sep = ifelse(args$sep == "tab", "\t", args$sep),
  check.names = FALSE,
  row.names = 1
)

pdf(file = args$output_pdf, paper = "a4")
if (args$package == "scater") {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(data.counts)))
  if (args$is_mito == TRUE) {
    ##retrieve mitochondrial genes in the dataset
    mito.genes <- grep(args$header_mito, rownames(sce))
    if (length(mito.genes) == 0)
      stop(paste(
        "No genes match mitochondrial pattern :",
        args$header_mito
      ))
    ##calculate QC metrics
    sce <-
     scater::calculateQCMetrics(sce,
                                feature_controls = list(Mito = mito.genes))
    ##search which cells had a high level of expressed mitochondrial genes
    high.mito <-
      scater::isOutlier(SingleCellExperiment::colData(sce)[, "pct_counts_Mito"],
                                   nmads = 3,
                                   type = "higher")
    ##remove those cells
    sce <- sce[, !high.mito]
    ##some QC plots
    hist(
      SingleCellExperiment::colData(sce)[, "total_counts"],
      breaks = 20,
      col = "grey80",
      main = "",
      xlab = "Log-total UMI count"
    )
    hist(
      SingleCellExperiment::colData(sce)[, "log10_total_features"],
      breaks = 20,
      col = "grey80",
      main = "",
      xlab = "Log-total number of expressed features"
    )
    hist(
      SingleCellExperiment::colData(sce)[, "pct_counts_Mito"],
      breaks = 20,
      col = "grey80",
      main = "",
      xlab = "Proportion of counts in mitochondrial genes"
    )
    ##Inspecting the most highly expressed genes
    print(scater::plotQC(sce,
                         type = "highest-expression",
                         n = ifelse(nrow(sce) < 50, nrow(sce), 50)))
  } else {
    ##calculate QC metrics
    sce <- scater::calculateQCMetrics(sce)
    ##search which cells had a low counts
    libsize.drop <-
      scater::isOutlier(
        SingleCellExperiment::colData(sce)[, "total_counts"],
        nmads = 3,
        type = "lower",
        log = TRUE
      )
    ##search for cells with few detected genes
    feature.drop <-
      scater::isOutlier(
        SingleCellExperiment::colData(sce)[, "total_features"],
        nmads = 3,
        type = "lower",
        log = TRUE
      )
    ##remove low quality cells
    sce <- sce[, !(libsize.drop | feature.drop)]
    cat("Remaining", ncol(sce), "cells.")
    ##Some QC plots
    hist(
      SingleCellExperiment::colData(sce)[, "total_counts"],
      breaks = 20,
      col = "grey80",
      main = "",
      xlab = "Log-total UMI count"
    )
    hist(
      SingleCellExperiment::colData(sce)[, "log10_total_features"],
      breaks = 20,
      col = "grey80",
      main = "",
      xlab = "Log-total number of expressed features"
    )
    ##Verify that the frequency of expression
    ##and the mean are positively correlated
    print(scater::plotQC(sce, type = "exprs-freq-vs-mean"))
    ##Inspecting the most highly expressed genes
    print(scater::plotQC(sce,
                         type = "highest-expression",
                         n = ifelse(nrow(sce) < 50, nrow(sce), 50)))
  }
  if (args$gene_filter == "min.cells") {
    numcells <- scater::nexprs(sce, byrow = TRUE)
    ##Filter genes detected in less than n cells
    numcells2 <- numcells >= args$min_cells
    sce <- sce[numcells2, ]
    cat("Keep", nrow(sce), "genes.")
  }
  if (args$gene_filter == "low.abundances") {
    ave.counts <- scater::calcAverage(sce)
    num.cells <- scater::nexprs(sce, byrow = TRUE)
    smoothScatter(
      log10(ave.counts),
      num.cells,
      ylab = "Number of cells",
      xlab = expression(Log[10] ~ "average count")
    )
    hist(
      log10(ave.counts),
      breaks = 20,
      main = "",
      col = "grey80",
      xlab = expression(Log[10] ~ "average count")
    )
    abline(
      v = log10(args$min_cells),
      col = "red",
      lwd = 2,
      lty = 2
    )
    to.keep <- ave.counts >= args$min_cells
    sce <- sce[to.keep, ]
    cat("Keep", nrow(sce), "genes.")
  }
}
if (args$package == "Seurat") {
  sce <- Seurat::CreateSeuratObject(
    raw.data =  data.counts,
    min.cells = args$min_cells,
    min.genes = args$min_genes
  )
  if (args$is_mito == TRUE) {
    ##retrieve mitochondrial genes in the dataset
    mito.genes <- grep(args$header_mito, rownames(sce@raw.data))
    if (length(mito.genes) == 0)
      stop(paste(
        "No genes match mitochondrial pattern :",
        args$header_mito
      ))
    ##calculate QC metrics
    percent.mito <- Matrix::colSums(sce@raw.data[mito.genes, ]) /
       Matrix::colSums(sce@raw.data)
    sce <- Seurat::AddMetaData(object = sce,
                               metadata = percent.mito,
                               col.name = "percent.mito")
    ##Filter low quality cells
    sce <- Seurat::FilterCells(sce,
                               subset.names = "percent.mito",
                               high.thresholds = 0.2)
    ##QC plot after filtering
    print(Seurat::VlnPlot(
      object = sce,
      features.plot = c("nGene", "nUMI", "percent.mito"),
      nCol = 3
    ))
  } else {
    print(Seurat::VlnPlot(
      object = sce,
      features.plot = c("nGene", "nUMI"),
      nCol = 2
    ))
  }
  ##Some QC plots
  hist(
    sce@meta.data[, "nUMI"],
    breaks = 20,
    main = "",
    col = "grey80",
    xlab = "Number of UMI/Cells"
  )
  hist(
    sce@meta.data[, "nGene"],
    breaks = 20,
    col = "grey80",
    main = "",
    xlab = "Number of genes detected/Cells"
  )
}
dev.off()

if (args$package == "scater") {
  filtered_data_counts <- sce@assays$data$counts
} else filtered_data_counts <- as.matrix(sce@data)

write.table(
  filtered_data_counts,
  file = args$output_matrix,
  sep = "\t",
  quote = F,
  col.names = args$colnames,
  row.names = T
)
