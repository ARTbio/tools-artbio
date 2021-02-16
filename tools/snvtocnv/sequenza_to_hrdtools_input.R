' sequenza_to_hrdtools_input.R
Takes data in the format of Sequenza output segments and
reformats it to be suitable for HRDtools input.
The output segments file contains one segment per row,
with the columns chr, start, end, copynumber, lohtype.
Usage: sequenza_to_hrdtools_input.R -i INPUT -s SOLUTIONS -o OUTPUT
Options:
  -i --input INPUT          Path to Sequenza output segments file
  -s --solutions SOLUTIONS  Path to Sequenza outputted list of alternative solutions
  -o --output OUTPUT        HRDtools input file
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

sequenza_data <- read_tsv(args[['--input']])
solutions_data <- read_tsv(args[['--solutions']])

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
  write_tsv(args[['--output']])

message(sprintf('Output written to %s', args[['output']]))

