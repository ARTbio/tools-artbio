# load packages that are provided in the conda env
options(show.error.messages = F,
       error = function() {
           cat(geterrmessage(), file = stderr()); q("no", 1, F)})
Sys.setenv(TZ = "Pacific/Auckland") # turnaround the tidyverse bug "In OlsonNames() : no Olson database found"


library(optparse)
library(sequenza)
library(BiocParallel)
library(tidyverse)

option_list <- list(
  make_option(
    c("-i", "--input"),
    default = NA,
    type = "character",
    help = "Path to Sequenza seqz processed segments file"
  ),
  make_option(
    c("-O", "--output_dir"),
    default = NA,
    type = "character",
    help = "Output directory"
  ),
  make_option(
    c("-s", "--sample_name"),
    default = NA,
    type = "character",
    help = "Sample name"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

data_file <- opt$input
output_dir <- opt$output_dir
sample_name <- opt$sample_name


## Processing seqz files : normalisation and segmentation for chromosomes 1 to 22
message(sprintf("\nExtraction step for %s", data_file))

chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
               "chr18", "chr19", "chr20", "chr21", "chr22")

segfile <- sequenza.extract(data_file,
                            verbose = TRUE,
                            chromosome.list = chrom_list)

## Estimation of cellularity and ploidy

segfile_cp <- sequenza.fit(segfile)
message(sprintf("\nEstimation step for %s\n", data_file))

## writing files and plots using default parameters
message(sprintf("\nWriting files and plots for %s\n", data_file))

sequenza.results(sequenza.extract = segfile,
                 cp.table = segfile_cp,
                 sample.id = sample_name,
                 out.dir = output_dir)
message(sprintf("\nOutput written to %s\n", output_dir))
