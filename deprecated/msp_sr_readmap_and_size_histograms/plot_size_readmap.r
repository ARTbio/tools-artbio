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
    make_option(c("-r", "--readmap_tab"), type="character", help="Path to file with tabular readmap"),
    make_option(c("-s", "--size_distribution_tab"), type="character", help="Path to file with tabular size distribution"),
    make_option("--readmap_pdf", type="character", help="Path to file with readmap plot"),
    make_option("--size_distribution_pdf", type="character", help="Path to file with size distribution plot"),
    make_option("--combi_pdf", type="character", help="Path to file with size distribution and readmap plot"),
    make_option("--title", type="character", help="Title for readmaps and size distribution"),
    make_option("--xlabel", type="character", help="xlabel for readmaps and size distribution"),
    make_option("--ylabel", type="character", help="ylabel for readmaps and size distribution"),
    make_option("--yrange", type="integer", help="Y-axis range"),
    make_option("--rows_per_page", type="integer", help="rows_per_page")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

## data frames implementation

rm=read.delim(args$readmap_tab, header=T, row.names=NULL)
n_samples=length(unique(rm$sample))
genes=unique(levels(rm$gene))
per_gene_readmap=lapply(genes, function(x) subset(rm, gene==x))
n_genes=length(per_gene_readmap)

size=read.delim(args$size_distribution_tab, header=T, row.names=NULL)
per_gene_size=lapply(genes, function(x) subset(size, gene==x))

## end of data frames implementation

## functions

plot_readmap=function(df, ...) {
combineLimits(xyplot(count~coord|factor(sample, levels=unique(sample))+reorder(gene, count, function(x) -sum(abs(x))),
data=df,
type='h',
scales= list(relation="free", x=list(rot=0, cex=0.7, axs="i", tck=0.5), y=list(tick.number=4, rot=90, cex=0.7)),
xlab=NULL, main=NULL, ylab=NULL,
as.table=T,
origin = 0,
horizontal=FALSE,
group=polarity,
col=c("red","blue"),
par.strip.text = list(cex=0.7),
...))
}

plot_size_distribution= function(df, ...) {
  smR.prepanel=function(x,y,...){; yscale=c(-max(abs(y)), max(abs(y)));list(ylim=yscale);}
  bc= barchart(count~as.factor(size)|factor(sample, levels=unique(sample))+gene, data = df, origin = 0,
    horizontal=FALSE,
group=polarity,
stack=TRUE,
    col=c('red', 'blue'),
    cex=0.75,
    scales=list(y=list(tick.number=4, rot=90, relation="free", cex=0.7), x=list(cex=0.7) ),
    prepanel=smR.prepanel,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    as.table=TRUE,
    newpage = T,
    par.strip.text = list(cex=0.7),
    ...)
  combineLimits(bc)
  }

## end of functions

## function parameters'

par.settings.readmap=list(layout.heights=list(top.padding=0, bottom.padding=-2.5), strip.background = list(col=c("lightblue","lightgreen")) )
par.settings.size=list(layout.heights=list(top.padding=-1, bottom.padding=-2.5), strip.background = list(col=c("lightblue","lightgreen")) )
par.settings.combination.readmap=list(layout.heights=list(top.padding=0, bottom.padding=-3), strip.background=list(col=c("lightblue","lightgreen")) )
par.settings.combination.size=list(layout.heights=list(top.padding=-2, bottom.padding=-0.5), strip.background=list(col=c("lightblue", "lightgreen")) )

## end of function parameters'

## GRAPHS

if (n_genes > 7) {page_height_simple = 11.69; page_height_combi=11.69; rows_per_page=args$rows_per_page} else {
                 rows_per_page= n_genes; page_height_simple = 2.5*n_genes; page_height_combi=page_height_simple*2 }
if (n_samples > 4) {page_width = 8.2677*n_samples/4} else {page_width = 8.2677*n_samples/3} # to test


pdf(file=args$readmap_pdf, paper="special", height=page_height_simple, width=page_width)
for (i in seq(1,n_genes,rows_per_page)) {
start=i
end=i+rows_per_page-1
if (end>n_genes) {end=n_genes}
if (args$yrange == 0) { readmap_plot.list=lapply(per_gene_readmap[start:end], function(x) plot_readmap(x, par.settings=par.settings.readmap)) } else {
readmap_plot.list=lapply(per_gene_readmap[start:end], function(x) plot_readmap(x, ylim=c(-args$yrange, args$yrange) , par.settings=par.settings.readmap)) }
args_list=c(readmap_plot.list, list(nrow=rows_per_page, ncol=1,
                                    top=textGrob("Read Maps (nucleotide coordinates)", gp=gpar(cex=1), just="top"),
                                    left=textGrob(args$ylabel, gp=gpar(cex=1), vjust=1, rot=90)
                                    )
           )
do.call(grid.arrange, args_list)
}
devname=dev.off()

pdf(file=args$size_distribution_pdf, paper="special", height=page_height_simple, width=page_width)
for (i in seq(1,n_genes,rows_per_page)) {
start=i
end=i+rows_per_page-1
if (end>n_genes) {end=n_genes}
plot.list=lapply(per_gene_size[start:end], function(x) plot_size_distribution(x, par.settings=par.settings.size) )
args_list=c(plot.list, list(nrow=rows_per_page, ncol=1,
                            top=textGrob("Size distributions (in nucleotides)", gp=gpar(cex=1), just="top"),
                            left=textGrob(args$ylabel, gp=gpar(cex=1), vjust=1, rot=90)
                            )
            )
do.call(grid.arrange, args_list)
}
devname=dev.off()

pdf(file=args$combi_pdf, paper="special", height=page_height_combi, width=page_width)
if (rows_per_page %% 2 != 0) { rows_per_page = rows_per_page + 1}
for (i in seq(1,n_genes,rows_per_page/2)) {
start=i
end=i+rows_per_page/2-1
if (end>n_genes) {end=n_genes}
if (args$yrange == 0) {readmap_plot.list=lapply(per_gene_readmap[start:end], function(x) plot_readmap(x, par.settings=par.settings.readmap)) } else {
readmap_plot.list=lapply(per_gene_readmap[start:end], function(x) plot_readmap(x, ylim=c(-args$yrange, args$yrange), par.settings=par.settings.readmap)) }
size_plot.list=lapply(per_gene_size[start:end], function(x) plot_size_distribution(x, strip=FALSE, par.settings=par.settings.combination.size))
plot.list=rbind(readmap_plot.list, size_plot.list )
args_list=c(plot.list, list(nrow=rows_per_page+1, ncol=1,
                            top=textGrob(args$title, gp=gpar(cex=1), just="top"),
                            left=textGrob(args$ylabel, gp=gpar(cex=1), vjust=1, rot=90),
                            sub=textGrob(args$xlabel, gp=gpar(cex=1), just="bottom")
                            )
            )
do.call(grid.arrange, args_list)
}
devname=dev.off()
