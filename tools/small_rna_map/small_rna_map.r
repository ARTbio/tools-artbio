library(optparse)

option_list <- list(
    make_option(c("-r", "--output_tab"), type="character", help="path to tabular file"),
    make_option("--output_pdf", type = "character", help="path to the pdf file with plot")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)

Table = read.delim(args$output_tab, header=T, row.names=NULL)

pdf(file=args$output_pdf, paper="special", height= 11.69, width = 8.2677)
plot(Table)
dev.off()
