#!/usr/bin/env Rscript

# A command-line interface to edgeR for use with Galaxy edger-repenrich
# written by Christophe Antoniewski drosofff@gmail.com 2017.05.30


# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# To not crash galaxy with an UTF8 error with not-US LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
library("tools")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "quiet", "q", 0, "logical",
  "help", "h", 0, "logical",
  "outfile", "o", 1, "character",
  "countsfile", "n", 1, "character",
  "factorName", "N", 1, "character",
  "levelNameA", "A", 1, "character",
  "levelNameB", "B", 1, "character",
  "levelAfiles", "a", 1, "character",
  "levelBfiles", "b", 1, "character",
  "alignmentA", "i", 1, "character",
  "alignmentB", "j", 1, "character",
  "plots" , "p", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# enforce the following required arguments
if (is.null(opt$outfile)) {
  cat("'outfile' is required\n")
  q(status=1)
}
if (is.null(opt$levelAfiles) | is.null(opt$levelBfiles)) {
  cat("input count files are required for both levels\n")
  q(status=1)
}
if (is.null(opt$alignmentA) | is.null(opt$alignmentB)) {
  cat("total aligned read files are required for both levels\n")
  q(status=1)
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
} else {
  FALSE
}

suppressPackageStartupMessages({
  library("edgeR")
  library("limma")
})

# build levels A and B file lists

library("rjson")
filesA <- fromJSON(opt$levelAfiles, method = "C", unexpected.escape = "error")
filesB <- fromJSON(opt$levelBfiles, method = "C", unexpected.escape = "error")
listA <- list()
indice = 0
listA[["level"]] <- opt$levelNameA
for (file in filesA) {
    indice = indice +1
    listA[[paste0(opt$levelNameA,"_",indice)]] <- read.delim(file, header=FALSE)
    }
listB <- list()
indice = 0
listB[["level"]] <- opt$levelNameB
for (file in filesB) {
    indice = indice +1
    listB[[paste0(opt$levelNameB,"_",indice)]] <- read.delim(file, header=FALSE)
    }

# build a counts table
counts <- data.frame(row.names=listA[[2]][,1])
for (element in names(listA[-1])) {
    counts<-cbind(counts, listA[[element]][,4])
    } 
for (element in names(listB[-1])) {
    counts<-cbind(counts, listB[[element]][,4])
    }
colnames(counts)=c(names(listA[-1]), names(listB[-1]))

# build aligned counts vector

filesi <- fromJSON(opt$alignmentA, method = "C", unexpected.escape = "error")
filesj <- fromJSON(opt$alignmentB, method = "C", unexpected.escape = "error")
sizes <- c()
for (file in filesi) {
    sizes <- c(sizes, read.delim(file, header=FALSE)[1,1])
    }
for (file in filesj) {
    sizes <- c(sizes, read.delim(file, header=FALSE)[1,1])
    }

# build a meta data object

meta <- data.frame(
    row.names=colnames(counts),
    condition=c(rep(opt$levelNameA,length(filesA)), rep(opt$levelNameB,length(filesB)) ),
    libsize=sizes
)


# Define the library size and conditions for the GLM
libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(meta$condition)


# Build a DGE object for the GLM
y <- DGEList(counts=counts, lib.size=libsize)

# Normalize the data
y <- calcNormFactors(y)
y$samples
# plotMDS(y) latter

# Estimate the variance
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
# plotBCV(y) latter

# Builds and outputs an object to contain the normalized read abundance in counts per million of reads
cpm <- cpm(y, log=FALSE, lib.size=libsize)
cpm <- as.data.frame(cpm)
colnames(cpm) <- colnames(counts)
if (!is.null(opt$countsfile)) {
    normalizedAbundance <- data.frame(Tag=rownames(cpm))
    normalizedAbundance <- cbind(normalizedAbundance, cpm)
    write.table(normalizedAbundance, file=opt$countsfile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

# Conduct fitting of the GLM
yfit <- glmFit(y, design)

# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

# Make the comparisons for the GLM
my.contrasts <- makeContrasts(
    paste0(opt$levelNameA,"_",opt$levelNameB," = ", opt$levelNameA, " - ", opt$levelNameB),
    levels = design
)

# Define the contrasts used in the comparisons
allcontrasts =  paste0(opt$levelNameA," vs ",opt$levelNameB)

# Conduct a for loop that will do the fitting of the GLM for each comparison
# Put the results into the results objects
    lrt <- glmLRT(yfit, contrast=my.contrasts[,1])
    plotSmear(lrt, de.tags=rownames(y))
    title(allcontrasts)
    res <- topTags(lrt,n=dim(c)[1],sort.by="none")$table
    results <- cbind(results,res[,c(1,5)])
    logfc <- cbind(logfc,res[c(1)])

# Add the repeat types back into the results.
# We should still have the same order as the input data
results$class <- listA[[2]][,2]
results$type <- listA[[2]][,3]

# Sort the results table by the FDR
results <- results[with(results, order(FDR)), ]

# Save the results
write.table(results, opt$outfile, quote=FALSE, sep="\t", col.names=FALSE)

# Plot Fold Changes for repeat classes and types

# open the device and plots
if (!is.null(opt$plots)) {
    if (verbose) cat("creating plots\n")
    pdf(opt$plots)
    plotMDS(y, main="Multidimensional Scaling Plot Of Distances Between Samples")
    plotBCV(y, xlab="Gene abundance (Average log CPM)", main="Biological Coefficient of Variation Plot")
    logFC <- results[, "logFC"]
    # Plot the repeat classes
    classes <- with(results, reorder(class, -logFC, median))
    par(mar=c(6,10,4,1))
    boxplot(logFC ~ classes, data=results, outline=FALSE, horizontal=TRUE,
        las=2, xlab="log(Fold Change)", main=paste0(allcontrasts, ", by Class"))
    abline(v=0)
    # Plot the repeat types
    types <- with(results, reorder(type, -logFC, median))
    boxplot(logFC ~ types, data=results, outline=FALSE, horizontal=TRUE,
        las=2, xlab="log(Fold Change)", main=paste0(allcontrasts, ", by Type"))
    abline(v=0)
}

# close the plot device
if (!is.null(opt$plots)) {
  cat("closing plot device\n")
  dev.off()
}

cat("Session information:\n\n")

sessionInfo()

