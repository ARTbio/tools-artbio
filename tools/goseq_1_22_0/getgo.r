suppressWarnings(suppressMessages(library("goseq")))
suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library("rtracklayer")))
suppressWarnings(suppressMessages(library("reshape2")))
sink(stdout(), type = "message")

option_list <- list(
    make_option(c("-gtf", "--gtf"), type="character", help = "Path to GTF file for which to fetch GO data"),
    make_option(c("-g", "--genome"), type="character", help = "Genome [used for looking up GO categories]"),
    make_option(c("-i", "--gene_id"), type="character", help="Gene ID format"),
    make_option(c("-c", "--cats"), type="character", help="Comma-seperated list of categories to fetch"),
    make_option(c("-o", "--output"), type="character", help="Path to output file")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:

gtf = args$gtf
genome = args$genome
gene_id = args$gene_id
output = args$output
cats = unlist(strsplit(args$cats, ','))
genes = unique(import.gff(gtf)$gene_id)
go_categories = getgo(genes, genome, gene_id, fetch.cats=cats)

# transform go category list to sth. more manipulatable in galaxy
go_categories <- lapply(go_categories, unlist)
go_categories = goseq:::reversemapping(go_categories)
go_categories = melt(go_categories)
colnames(go_categories) = c("#gene_id", "go_category")

write.table(go_categories, output, sep="\t", row.names = FALSE, quote = FALSE)
sessionInfo()