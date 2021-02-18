# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

library(optparse)
library(tidyverse)

option_list = list(
  make_option(
    c("-i", "--input"),
    default = NA,
    type = "character",
    help = "Path to Sequenza output segments file"
  ),
  make_option(
    c("-o", "--output"),
    default = NA,
    type = "character",
    help = "output file, to be used as input for HRDetect"
  ),
  make_option(
    c("-s", "--solutions"),
    default = NA,
    type = "character",
    help = "Path to Sequenza list of alternative solutions"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))


sequenza_data <- read_tsv(opt$input)
solutions_data <- read_tsv(opt$solutions)

ploidy <- round(solutions_data$ploidy[1])
cellularity <- solutions_data$cellularity[1]

reformatted <- sequenza_data %>%
  select(
    chr=chromosome,
    start=start.pos,
    end=end.pos,
    copynumber=CNt,
    A, B
  ) %>%
  mutate(
    ploidy=ploidy,
    cellularity=cellularity,
    lohtype=case_when(
      copynumber==0 ~ 'HOMD',
      B==0 & A==ploidy ~ 'NLOH',
      B==0 & A<ploidy & A > 0 ~ 'DLOH',
      copynumber>ploidy & A>B ~ 'ASCNA',
      copynumber>ploidy & A==B ~ 'BCNA',
      TRUE ~ 'HET'
    )
  )

message('Preview of output:')
print(reformatted)

reformatted %>%
  write_tsv(opt$output)

message(sprintf('Output written to %s', opt$output))

