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

for(i in 1:nrow(Table)){
  if (Table$Polarity[i] == "R"){
	Table$Nbr_reads[i]= -1 * Table$Nbr_reads[i]
	}
}

i=1
plots <- list()
for (read in unique(Table$Chromosome,Table$Dataset)) {
  
  bp <- ggplot(Table[Table$Chromosome==read,], 
  aes(x=Coordinate, y=Nbr_reads, fill=Polarity))+ 
  geom_bar(stat="identity", position=position_dodge(),width=1, color= "black" )+
  facet_wrap(Dataset~Chromosome)+
  geom_hline(yintercept=0)+
  geom_text(aes(label=Nbr_reads), position = position_dodge(width = 0.9),vjust=1.5, size=2) 
  plots[[i]] <- bp+scale_x_continuous(limits=c(0,Table$Chrom_length[i]))
  i=i+1
}

p<-do.call(marrangeGrob,list(grobs=plots,ncol=1,nrow=3));
ggsave(file=args$output_pdf, plot=p, height=11.69, width=8.2677)



