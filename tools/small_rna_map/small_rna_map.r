library(optparse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gtable)

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

# The second plot
p2 <- ggplot(Table, aes(x = Coordinate, y = Median)) +
  geom_point(aes(col=Median),size = 1, colour="black") + 
  geom_line()+
  expand_limits(y = seq(0,max(Table$Median),by=4)) +
  facet_wrap(Dataset~Chromosome, scales="free") +
  geom_segment(aes(y = Nbr_reads, 
			 x = 0,
			 yend=Nbr_reads,
			 xend=Chrom_length), alpha=0
		   )+
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        legend.position = "bottom")

# Transforme ggplot graphs on list of graphs
plot.list1 <- by(data     = Table,
                INDICES  = c(Table$Chromosome),
                simplify = TRUE,
                FUN      = function(x) {
                  p %+% x 
                })
 
plot.list2 <- by(data     = Table,
                INDICES  = c(Table$Chromosome),
                simplify = TRUE,
                FUN      = function(x) {
                  p2 %+% x 
                })

dual_axis <- function(v1,v2){
# Get the ggplot grobs
			g1 <- ggplot_gtable(ggplot_build(v1))
			g2 <- ggplot_gtable(ggplot_build(v2))
# Get the locations of the plot panels in g1.
			pp <- c(subset(g1$layout, grepl("panel", g1$layout$name), se = t:r))
# Overlap panels for second plot on those of the first plot
			g <- gtable_add_grob(g1, g2$grobs[grepl("panel", g1$layout$name)], 
      			pp$t, pp$l, pp$b, pp$l)
		}
 
plots <- list()
 
 

len = length(plot.list1)
for(i in 1:len ) {plots[[i]] <- dual_axis(plot.list1[[i]],plot.list2[[i]])}
# Plotting in multiple pages with different rows
multi.plot<-do.call(marrangeGrob,list(grobs=plots,ncol=1,nrow=8,top=NULL, bottom="Coordinates(nt)", left="Number of reads", right = mylegend));
ggsave(args$output_pdf, device="pdf", plot=multi.plot, height=11.69, width=8.2)

