library(optparse)
library(ggplot2)
library(gridExtra)

option_list <- list(
    make_option(c("-r", "--output_tab"), type="character", help="path to tabular file"),
    make_option("--output_pdf", type = "character", help="path to the pdf file with plot")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)

theme_set(theme_bw())
Table = read.delim(args$output_tab, header=T, row.names=NULL)

Table <- within(Table, Nbr_reads[Polarity=="R"] <- (Nbr_reads[Polarity=="R"]*-1))


i=1
plots <- list()
for (read in unique(Table$Chromosome,Table$Dataset)) {
  
  bp <- ggplot(Table[Table$Chromosome==read,], 
  aes(x=Coordinate, y=Nbr_reads, color=Polarity))+ 
  geom_segment(aes(y = 0, 
                   x = Coordinate, 
                   yend = Nbr_reads, 
                   xend = Coordinate), size=1)+
  facet_wrap(Dataset~Chromosome, scales="free")+
  geom_hline(yintercept=0)
  
  plots[[i]] <- bp+scale_x_continuous(limits=c(0,Table$Chrom_length[i]))
  i=i+1
}

p<-do.call(marrangeGrob,list(grobs=plots,ncol=1,nrow=4,top=NULL));
ggsave(file=args$output_pdf, plot=p, height=11.69, width=8)



