library(optparse)
library(ggplot2)
library(gridExtra)

option_list <- list(
    make_option(c("-r", "--output_tab"), type="character", help="path to tabular file"),
    make_option("--output_pdf", type = "character", help="path to the pdf file with plot")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)

Table = read.delim(args$output_tab, header=T, row.names=NULL)


i=1
plots <- list()
for (read in unique(Table$Chromosome)) {
  bp <- ggplot(Table[Table$Chromosome==read,], 
  aes(x=Coordinate, y=Nbr_reads, fill=Polarity))+ 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=1)+
  facet_grid(Chromosome ~ Dataset)+
  theme(strip.text.x = element_text(size=8),
  strip.text.y = element_text(size=12, face="bold"),
  strip.background = element_rect(colour="white", fill="#CCCCFF"))+
  geom_text(aes(label=Nbr_reads), vjust=1.2, color="white", position =position_dodge(0.9),size=2)
  plots[[i]] <- bp+xlim(0,Table$Chrom_length[i])
  i=i+1
}
pdf(args$output_pdf, paper="special", height=11.69, width=8.2677)
do.call(marrangeGrob,list(grobs=plots,ncol=1,nrow=3));
dev.off()
