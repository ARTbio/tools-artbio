#!/usr/bin/env Rscript

# Help to be printed
hlp_description = "This script takes a dataframe containing the coverage values for different miRNAs and plots them."
hlp_usage = "Usage : coverage_ploting.R --dataframe [FILE] --type ['relative' or 'absolute'] --output [FILE]"
hlp_dataframe = "--dataframe\tFILE\tThis is a dataframe containing coverage values obtained from mircounts.py"
hlp_type = "--type\t\tSTRING\tType of ploting, either relative or absoute coverage values (default='relative')"
hlp_output = "--output\tFILE\tFile to output the pdf to\n"

hlp = paste(hlp_description,hlp_usage,hlp_dataframe,hlp_type,hlp_output,sep="\n")
#print(hlp)

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

library(optparse)
library(lattice)

# Get arguments
option_list <- list(
    make_option(c("-d", "--dataframe"), type="character",
                help="Dataframe containing coverage values obtained from mircounts.py"),
    make_option(c("-t", "--type"), type="character", default="relative",
                help="Type of plotting, either relative or absoute coverage values (default='relative')"),
    make_option(c("-o", "--output"), type="character", help="File to output the pdf to")
    )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

if ( !('dataframe' %in% names(args)) || !('output' %in% names(args))) {
    stop("'--dataframe' and '--output' parametters are not optional. Please retry.")
}

# Plot
coverage = read.delim(args$dataframe, header=T)
if (args$type =="relative") {
    graph = xyplot(Norm_count~Norm_offset | Mir_hairpin, data=coverage, col=c("darkblue"), type="l", lwd=1.5,
                   scales=list(x=list(cex=.5), y=list(cex=.5)), par.strip.text=list(cex=.5),
                   strip=strip.custom(which.given=1, bg="lightblue"), layout=c(4,15),
                   as.table=T, xlab="Normalized Counts", ylab="Normalized coordinates",
                   main="miRNA coverage maps")
} else {
    graph = xyplot(Count~Offset | Mir_hairpin, data=coverage, col=c("darkblue"), type="l", lwd=1.5,
                   scales=list(x=list(cex=.5), y=list(cex=.5)), par.strip.text=list(cex=.5),
                   strip=strip.custom(which.given=1, bg="lightblue"), layout=c(4,15),
                   as.table=T, xlab="Counts", ylab="Coordinates",
                   main="miRNA coverage plots")
}

# PDF output
pdf(file=args$output, paper="special", height=11.69, width=8.2677)
plot(graph, newpage = T)
dev.off()
