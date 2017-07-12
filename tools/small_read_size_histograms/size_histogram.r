## Setup R error handling to go to stderr
options( show.error.messages=F,
         error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(optparse)

# Parse arguments
option_list <- list(
    make_option(c("-g", "--global"), type="character", help="Whether distribution is plotted globally or by chromosome"),
    make_option(c("-s", "--size_distribution_tab"), type="character", help="Path to file with tabular size distribution"),
    make_option("--size_distribution_pdf", type="character", help="Path to file with size distribution plot"),
    make_option("--title", type="character", help="Title for readmaps and size distribution"),
    make_option("--ylabel", type="character", help="ylabel for readmaps and size distribution"),
    make_option("--yrange", type="integer", help="Y-axis range"),
    make_option("--rows_per_page", type="integer", help="rows_per_page")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

##cheetahtemplate data frame implementation
size=read.delim(args$size_distribution_tab, header=T, row.names=NULL)
n_samples = length(unique (size$sample))
n_genes = length (unique (levels(size$gene)))

if (args$yrange != 0) {
   # This is used for specifying the y-axis limits
   ylim=c(-args$yrange, args$yrange)
} else { ylim="" }

par.settings.size=list(layout.heights=list(top.padding=1, bottom.padding=1),
                       strip.background = list(col = c("lightblue", "lightgreen"))
                       )

smR.prepanel=function(x,y,...){; yscale=c(-max(abs(y)), max(abs(y)));list(ylim=yscale);} # use if one want y axis in the middle of the plot

plot_size_distribution = function(df, ...) {
   bc= barchart(count~as.factor(size)|factor(sample, levels=unique(sample))+gene, data = df, origin = 0,
                horizontal=FALSE,
            group=polarity,
            stack=TRUE,
                col=c('red', 'blue'),
                cex=0.75,
                scales=list(y=list(tick.number=4, rot=90, relation="free", cex=0.5, alternating=T), x=list(cex=.6 ) ),
                xlab = "readsize in nucleotides",
                ylab = args$ylabel,
                main = args$title,
                par.strip.text = list(cex=0.75),
                as.table=TRUE,
                newpage = T,
                ...)

    combineLimits(update(useOuterStrips(bc,
                                        strip.left = strip.custom(par.strip.text = list(cex=0.5))
                                        ),
                  layout=c(n_samples,args$rows_per_page)),
                  margin.x=F, margin.y=1)
    }

# per_gene_size=lapply(genes, function(x) subset(size, gene==x)) # no object in this script

if (args$global == "no") {
width = 8.2677*n_samples/4
} else { width = 8.2677 }

options(warn=-1)
pdf(file=args$size_distribution_pdf, paper="special", height=11.69, width=width)

if (ylim == "" && args$global=="no") {
    plot_size_distribution(size, par.settings=par.settings.size)
   }
if (ylim != "" && args$global=="no") { plot_size_distribution(size, par.settings=par.settings.size, ylim=ylim)
   }
if (ylim == "" && args$global=="yes") {  bc= barchart(count~as.factor(size)|factor(sample, levels=unique(sample)),
        data = size, origin = 0,
        horizontal=FALSE,
        group=polarity,
        stack=TRUE,
        col=c('red', 'blue'),
        scales=list(y=list(tick.number=4, rot=90, relation="same"), cex=1),
        xlab = "readsize in nucleotides",
        ylab = args$ylabel,
        main = args$title, as.table=TRUE, newpage = T,
        aspect=0.5,
        strip = strip.custom(par.strip.text = list(cex = 1), which.given=1, bg="lightblue")
        )
   bc
   }
if (ylim != "" && args$global=="yes") {  bc= barchart(count~as.factor(size)|factor(sample, levels=unique(sample)),
        data = size, origin = 0,
        horizontal=FALSE,
        group=polarity,
        stack=TRUE,
        col=c('red', 'blue'),
        scales=list(y=list(tick.number=4, rot=90, relation="same"), cex=1),
        xlab = "readsize in nucleotides",
        ylab = args$ylabel,
        ylim = ylim,
        main = args$title, as.table=TRUE, newpage = T,
        aspect=0.5,
        strip = strip.custom(par.strip.text = list(cex = 1), which.given=1, bg="lightblue")
        )
   bc
   }

devname=dev.off()
