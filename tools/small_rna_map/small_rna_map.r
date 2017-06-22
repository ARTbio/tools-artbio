library(optparse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

option_list <- list(
    make_option(c("-r", "--output_tab"), type="character", help="path to tabular file"),
    make_option("--output_pdf", type = "character", help="path to the pdf file with plot")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)

theme_set(theme_bw()) #a theme with a white background
Table = read.delim(args$output_tab, header=T, row.names=NULL)
Table <- within(Table, Nbr_reads[Polarity=="R"] <- (Nbr_reads[Polarity=="R"]*-1))

#To assign colors to categorical variables in ggplot2 that have stable mapping
myColors <- brewer.pal(3,"Set1")
names(myColors) <- levels(Table$Polarity)
colScale <- scale_colour_manual(name = "Polarity",values = myColors)

#Make initial figures
p <- ggplot(Table, aes(x=Coordinate, y=Nbr_reads, colour=Polarity)) +
  colScale+
  geom_segment(aes(y = 0, 
                   x = Coordinate, 
                   yend = Nbr_reads, 
                   xend = Coordinate,
			 color=Polarity),
			 alpha=1
               ) +
  geom_segment(aes(y = Nbr_reads, 
			 x = 0,
			 yend=Nbr_reads,
			 xend=Chrom_length), alpha=0
		   )+
  facet_wrap(Dataset~Chromosome, scales="free", nrow=1, labeller = label_wrap_gen(multi_line = FALSE))+
  scale_y_continuous(breaks = function(x) round(pretty(seq(-(max(x) + 1), (max(x) + 1)))))+#to display only integer values on y axis
  geom_hline(yintercept=0, size=0.3)+
  theme(strip.text = element_text(size = 6, lineheight = 0.1), #specify strip size
        panel.grid.major = element_line(colour = "#ffffff"),#conceal major grid lines
        panel.grid.minor = element_line(colour = "#ffffff"),#conceal minor grid lines
        axis.title = element_blank(),# Conceal axis titles
        axis.text = element_text(size = 6))#modify the size of tick labels along axes
#extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(p);

p <- p+ theme(legend.position = "none")# Hide the repeated caption
#transforme ggplot graphs on list of graphs
plot.list <- by(data     = Table,
                INDICES  = c(Table$Chromosome),
                simplify = TRUE,
                FUN      = function(x) {
                  p %+% x 
                })
#Plotting in multiple pages with different rows
multi.plot <- marrangeGrob(grobs = plot.list,nrow  = 8, ncol = 1, top=NULL, bottom="Coordinates(nt)", left="Number of reads", right = mylegend);
ggsave(args$output_pdf, device="pdf", plot=multi.plot, height=11.69, width=8.2)

