' segmentation_sequenza.R
Takes data in the format of Sequenza output segz and
reformats it to calculate lohtypes.
Creates an output directory with different files and plots
Necessary files for next steps : *alternative_solutions.txt and *segments.txt
Usage: segmentation_sequenza.R -i INPUT -s SAMPLE_NAME -d OUTPUT_DIR
Options:
-i --input INPUT          Path to Sequenza seqz processed segments file
-s --sample-name SAMPLE_NAME
-d --output-dir OUTPUT_DIR  Output directory

' -> doc

library(docopt)
library(sequenza)
library(BiocParallel)
library(tidyverse)

args <- docopt(doc)

data.file <- args[['--input']]
output.dir <- args[['--output-dir']]
sample.name <- args[['--sample-name']]

## Processing seqz files : normalisation and segmentation for chromosomes 1 to 22
message(sprintf('\nExtraction step for %s', data.file))

segfile <- sequenza.extract(data.file, 
                            verbose = TRUE,
                            chromosome.list = as.character(c(1:22)))

## Estimation of cellularity and ploidy 

segfile.CP <- sequenza.fit(segfile)
message(sprintf('\nEstimation step for %s\n', data.file))

## writing files and plots using default parameters
message(sprintf('\nWriting files and plots for %s\n', data.file))

sequenza.results(sequenza.extract = segfile,
                 cp.table = segfile.CP, 
                 sample.id = sample.name,
                 out.dir = output.dir )


message(sprintf('\nOutput written to %s\n', output.dir))



