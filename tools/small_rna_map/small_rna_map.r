library(optparse)
library(ggplot2)

option_list <- list(
    make_option(c("-r", "--output_tab"), type="character", help="path to tabular file"),
    make_option("--output_pdf", type = "character", help="path to the pdf file with plot")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)

Table = read.delim(args$output_tab, header=T, row.names=NULL)

pdf(file=args$output_pdf, paper="special", height= 11.69, width = 8.2677)
for (read in unique(Table$Chromosome)) {
          print(ggplot(Table[Table$Chromosome==read,], aes(x=Coordinate, y=Nbr_reads, fill=Polarity))+ 
          geom_bar(stat="identity", color="black", position=position_dodge())+
	    facet_grid(Chromosome ~ Dataset)+
          theme(strip.text.x = element_text(size=8, angle=75),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="red", fill="#CCCCFF"))+
          geom_text(aes(label=Nbr_reads), vjust=1.2, color="white", position =position_dodge(0.9),size=3)) 
	}

dev.off()
