## Setup R error handling to go to stderr
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# options(warn = -1)
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(optparse)

option_list <- list(
    make_option(c("-f", "--first_dataframe"), type="character", help="path to first dataframe"),
    make_option(c("-e", "--extra_dataframe"), type="character", help="path to additional dataframe"),
    make_option("--first_plot_method", type = "character", help="How additional data should be plotted"),
    make_option("--extra_plot_method", type = "character", help="How additional data should be plotted"),
    make_option("--output_pdf", type = "character", help="path to the pdf file with plots")
    )
 
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)
 
# data frames implementation
## first table
Table = read.delim(args$first_dataframe, header=T, row.names=NULL)
if (args$first_plot_method == "Counts" | args$first_plot_method == "Size") {
    Table <- within(Table, Counts[Polarity=="R"] <- (Counts[Polarity=="R"]*-1))
}
n_samples=length(unique(Table$Dataset))
genes=unique(levels(Table$Chromosome))
per_gene_readmap=lapply(genes, function(x) subset(Table, Chromosome==x))
per_gene_limit=lapply(genes, function(x) c(1, unique(subset(Table, Chromosome==x)$Chrom_length)) )
n_genes=length(per_gene_readmap)
# second table
if (args$extra_plot_method != '') {
    ExtraTable=read.delim(args$extra_dataframe, header=T, row.names=NULL)
    if (args$extra_plot_method == "Counts" | args$extra_plot_method=='Size') {
        ExtraTable <- within(ExtraTable, Counts[Polarity=="R"] <- (Counts[Polarity=="R"]*-1))
        }
    per_gene_size=lapply(genes, function(x) subset(ExtraTable, Chromosome==x))
    }

## functions

plot_unit = function(df, method=args$first_plot_method, ...) {
    if (method == 'Counts') {
        p = xyplot(Counts~Coordinate|factor(Dataset, levels=unique(Dataset))+factor(Chromosome, levels=unique(Chromosome)),
        data=df,
        type='h',
        lwd=1.5,
        scales= list(relation="free", x=list(rot=0, cex=0.7, axs="i", tck=0.5), y=list(tick.number=4, rot=90, cex=0.7)),
        xlab=NULL, main=NULL, ylab=NULL,
        as.table=T,
        origin = 0,
        horizontal=FALSE,
        group=Polarity,
        col=c("red","blue"),
        par.strip.text = list(cex=0.7),
        ...)
    } else if (method != "Size") {
        p = xyplot(eval(as.name(method))~Coordinate|factor(Dataset, levels=unique(Dataset))+factor(Chromosome, levels=unique(Chromosome)),
        data=df,
        type='p',
        pch=19,
        cex=0.35,
        scales= list(relation="free", x=list(rot=0, cex=0.7, axs="i", tck=0.5), y=list(tick.number=4, rot=90, cex=0.7)),
        xlab=NULL, main=NULL, ylab=NULL,
        as.table=T,
        origin = 0,
        horizontal=FALSE,
        group=Polarity,
        col=c("red","blue"),
        par.strip.text = list(cex=0.7),
        ...)
    } else {
        p = barchart(Counts~as.factor(Size)|factor(Dataset, levels=unique(Dataset))+Chromosome, data = df, origin = 0,
                     horizontal=FALSE,
                     group=Polarity,
                     stack=TRUE,
                     col=c('red', 'blue'),
                     scales=list(y=list(tick.number=4, rot=90, relation="free", cex=0.7), x=list(rot=0, cex=0.7, axs="i", tck=0.5)),
        xlab = NULL,
        ylab = NULL,
        main = NULL,
        as.table=TRUE,
        newpage = T,
        par.strip.text = list(cex=0.7),
        ...)
    }
    if (args$extra_plot_method != '') {
    combineLimits(p)
    }
}

plot_single <- function(df, method=args$first_plot_method, ...) {
        if (method == 'Counts') {
    
    } else {
        p= barchart(Counts~as.factor(Size)|factor(Dataset, levels=unique(Dataset))+Chromosome, data = df, origin = 0,
                    horizontal=FALSE,
                    group=Polarity,
                    stack=TRUE,
                    col=c('red', 'blue'),
                    scales=list(y=list(tick.number=4, rot=90, relation="free", cex=0.5, alternating=T), x=list(rot=0, cex=0.6, tck=0.5)),
                    xlab=NULL, main=NULL, ylab=NULL,
                    par.strip.text = list(cex=0.75),
                    as.table=TRUE,
                    ...)
          p = update(useOuterStrips(p, strip.left=strip.custom(par.strip.text = list(cex=0.5))))
          p = combineLimits(p, extend=TRUE, margin.x=1, margin.y=1)
          # p = update(p, layout=c(n_samples, 7))
          # p = update( useOuterStrips(combineLimits(p), strip.left=strip.custom(par.strip.text = list(cex=0.5))), layout=c(n_samples, 7) )

          # to avoid irritating message In heights.x[pos.heights[[nm]]] <- heights.settings[[nm]] * heights.defaults[[nm]]$x :
          # number of items to replace is not a multiple of replacement length
          # uncomment "options(warn = -1)" for effect
        }
}

## function parameters

par.settings.firstplot = list(layout.heights=list(top.padding=11, bottom.padding = -14))
par.settings.secondplot=list(layout.heights=list(top.padding=11, bottom.padding = -15), strip.background=list(col=c("lavender","deepskyblue")))
#par.settings.single_plot=list(strip.background=list(col=c("lavender","deepskyblue"))) # layout.heights=list(top.padding=-11, bottom.padding=-11),
par.settings.single_plot=list(layout.heights=list(top.padding=-1, bottom.padding = -5),strip.background = list(col = c("lightblue", "lightgreen")) ) # layout.heights=list(top.padding=1, bottom.padding=1), 
title_first_method = list(Counts="Read Counts", Coverage="Coverage depths", Median="Median sizes", Mean="Mean sizes", Size="Size Distributions")
title_extra_method = list(Counts="Read Counts", Coverage="Coverage depths", Median="Median sizes", Mean="Mean sizes", Size="Size Distributions")
legend_first_method =list(Counts="Read count", Coverage="Coverage depth", Median="Median size", Mean="Mean size", Size="Read count")
legend_extra_method =list(Counts="Read count", Coverage="Coveragedepth", Median="Median size", Mean="Mean size", Size="Read count")
bottom_first_method =list(Counts="Coordinates (nbre of bases)",Coverage="Coordinates (nbre of bases)", Median="Coordinates (nbre of bases)", Mean="Coordinates (nbre of bases)", Size="Sizes of reads")
bottom_extra_method =list(Counts="Coordinates (nbre of bases)",Coverage="Coordinates (nbre of bases)", Median="Coordinates (nbre of bases)", Mean="Coordinates (nbre of bases)", Size="Sizes of reads")

## Plotting Functions

double_plot <- function(...) {
    if (n_genes > 5) {page_height=15; rows_per_page=10} else {
                     rows_per_page= 2 * n_genes; page_height=1.5*n_genes}
    if (n_samples > 4) {page_width = 8.2677*n_samples/4} else {page_width = 7 * n_samples/2}
    pdf(file=args$output_pdf, paper="special", height=page_height, width=page_width)
    for (i in seq(1,n_genes,rows_per_page/2)) {
        start=i
        end=i+rows_per_page/2-1
        if (end>n_genes) {end=n_genes}
        first_plot.list = lapply(per_gene_readmap[start:end], function(x) plot_unit(x, strip=FALSE, par.settings=par.settings.firstplot))
        second_plot.list = lapply(per_gene_size[start:end], function(x) plot_unit(x, method=args$extra_plot_method, par.settings=par.settings.secondplot))
        plot.list=rbind(second_plot.list, first_plot.list)
        args_list=c(plot.list, list(nrow=rows_per_page+1, ncol=1,
                                    top=textGrob(paste(title_first_method[[args$first_plot_method]], "and", title_extra_method[[args$extra_plot_method]]), gp=gpar(cex=1), vjust=2, just="top"),
                                    left=textGrob(paste(legend_first_method[[args$first_plot_method]], "/", legend_extra_method[[args$extra_plot_method]]), gp=gpar(cex=1), vjust=2, rot=90),
                                    sub=textGrob(paste(bottom_first_method[[args$first_plot_method]], "/", bottom_extra_method[[args$extra_plot_method]]), gp=gpar(cex=1), just="bottom", vjust=4)
                                    )
                   )
        do.call(grid.arrange, args_list)
        }
    devname=dev.off()
}

single_plot <- function(...) {
    width = 8.2677 * n_samples / 2
    rows_per_page=8
    pdf(file=args$output_pdf, paper="special", height=11.69, width=width) #
    for (i in seq(1,n_genes,rows_per_page)) {
        start=i
        end=i+rows_per_page-1
        if (end>n_genes) {end=n_genes}
        plot.list = lapply(per_gene_readmap[start:end], function(x) plot_single(x, strip=FALSE, par.settings=par.settings.single_plot)) # strip=FALSE,
        args_list=c(plot.list, list( ncol=1, #  nrow=1,
                                    top=textGrob(title_first_method[[args$first_plot_method]], gp=gpar(cex=1), vjust=2, just="top"),
                                    left=textGrob(legend_first_method[[args$first_plot_method]], gp=gpar(cex=1), vjust=2, rot=90),
                                    sub=textGrob(bottom_first_method[[args$first_plot_method]], gp=gpar(cex=1), just="bottom", vjust=4)
                                    )
                     )
       # do.call(grid.arrange, c(plot.list, list(ncol=1, heights=rep(1,5))))
       do.call(grid.arrange, args_list)
        }
    dev.off()
}

# main

if (args$extra_plot_method != '') { double_plot() }
if (args$first_plot_method == "Size" & args$extra_plot_method == '') {
    single_plot()
}









