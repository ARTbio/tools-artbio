options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("goseq")
    library("optparse")
    library("rtracklayer")
    library("reshape2")
})

sink(stdout(), type = "message")

option_list <- list(
    # make_option(c("-gtf", "--gtf"), type="character", help = "Path to GTF file for which to fetch GO data"),
    make_option(c("-g", "--genome"), type="character", help = "Genome [used for looking up GO categories]"),
    make_option(c("-i", "--gene_id"), type="character", help="Gene ID format"),
    make_option(c("-c", "--cats"), type="character", help="Comma-seperated list of categories to fetch"),
    make_option(c("-o", "--output"), type="character", help="Path to output file")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# vars

gtf = args$gtf
genome = args$genome
gene_id = args$gene_id
output = args$output
cats = unlist(strsplit(args$cats, ','))

# retrieve and transform data
genes = unique(import.gff(gtf)$gene_id)
go_categories = getgo(genes, genome, gene_id, fetch.cats=cats)
go_categories = goseq:::reversemapping(go_categories)
go_categories = melt(go_categories)

write.table(go_categories, output, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
sessionInfo()