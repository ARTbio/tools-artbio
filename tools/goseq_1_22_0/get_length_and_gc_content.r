# originally by Devon Ryan, https://www.biostars.org/p/84467/
sink(stdout(), type = "message")

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(optparse)
library(data.table)

option_list <- list(
    make_option(c("-g","--gtf"), type="character", help="Input GTF file with gene / exon information."),
    make_option(c("-f","--fasta"), type="character", default=FALSE, help="Fasta file that corresponds to the supplied GTF."),
    make_option(c("-o","--output"), type="character", default=FALSE, help="Output file with gene name, length and GC content.")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

GTFfile = args$gtf
FASTAfile = args$fasta
output_file = args$output

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome=NA, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
output <- setDT(data.frame(output), keep.rownames = TRUE)[]
colnames(output) <- c("#gene_id", "length", "GC")

write.table(output, file=output_file, row.names=FALSE, quote=FALSE, sep="\t")

sessionInfo()