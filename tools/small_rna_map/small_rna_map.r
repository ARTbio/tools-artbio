library(optparse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gtable)
library(grid)
 
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
        axis.text = element_text(size = 6),#modify the size of tick labels along axes
        legend.position = "none") # Hide the repeate caption

# Create legend
mylegend <- legendGrob(c("F", "R", "Median", "Mean"), pch=22,
                     gp=gpar(col = c("red","blue","black","yellow"), fill = c("red","blue","black","yellow")))
 
 
# The second plot
cols<- c("Median"="#000000", "Mean"="#fffa00")
p2 <- ggplot(Table, aes(x = Coordinate, group=1)) +
  geom_point(aes(y=Median, colour="Median"), alpha=1, size = 1) +
  geom_point(aes(y=Mean, colour="Mean"), alpha= 0.3, size = 1.2)+
  scale_colour_manual(name="", values=cols)+ 
  expand_limits(y = seq(0,max(Table$Median),by=5)) +
  facet_wrap(Dataset~Chromosome, scales="free", nrow=1, labeller = label_wrap_gen(multi_line = FALSE))+
  geom_segment(aes(y = Nbr_reads, 
			 x = 0,
			 yend=Nbr_reads,
			 xend=Chrom_length), alpha=0
		   )+
  scale_y_continuous(limits = c(0,max(Table$Median)), position = "right")+
  theme(strip.text = element_text(size = 6, lineheight = 0.1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        axis.text = element_text(size = 6),
        legend.position = "none"
        )

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
 
# A function to get the original tick mark length
plot_theme <- function(p) {
  plyr::defaults(p$theme, theme_get())
}
 
# ggplot contains many labels that are themselves complex grob; 
# usually a text grob surrounded by margins.
# When moving the grobs from, say, the left to the right of a plot,
# Make sure the margins and the justifications are swapped around.
# The function below does the swapping.
# Taken from the cowplot package:
# https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
 
hinvert_title_grob <- function(grob){
 
  # Swap the widths
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
 
  # Fix the justification
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}
 
dual_axis <- function(v1,v2){
# Get the ggplot grobs
	g1 <- ggplot_gtable(ggplot_build(v1))
	g2 <- ggplot_gtable(ggplot_build(v2))
# Get the locations of the plot panels in g1.
	pp <- c(subset(g1$layout, grepl("panel", g1$layout$name), se = t:r))
# Overlap panels for second plot on those of the first plot
	g <- gtable_add_grob(g1, g2$grobs[grepl("panel", g1$layout$name)], 
      	pp$t, pp$l, pp$b, pp$l)
# Get the y axis from g2 (axis line, tick marks, and tick mark labels)
	index <- which(g2$layout$name == "axis-r-1-1")  # Which grob.   
	yaxis <- g2$grobs[[index]]      # Extract the grob
	ticks <- yaxis$children[[2]]  # swap tick marks and tick mark labels
 
# Move the tick marks, Tick mark lengths can change. 
	tml <- plot_theme(p)$axis.ticks.length   # Tick mark length
	ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + tml
 
# Swap margins and fix justifications for the tick mark labels
ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
 
# Put ticks back into yaxis
yaxis$children[[2]] <- ticks
 
# Put the transformed yaxis on the right side of g1
g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, yaxis, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, 
   clip = "off")
 
}
 
plots <- list()
len = length(plot.list1)
for(i in 1:len ) {plots[[i]] <- dual_axis(plot.list1[[i]],plot.list2[[i]])}
# Plotting in multiple pages with different rows
multi.plot<-do.call(marrangeGrob,list(grobs=plots,ncol=1,nrow=8,top=NULL, bottom="Coordinates(nt)", left="Number of reads", right= mylegend))
ggsave(args$output_pdf, device="pdf", plot=multi.plot, height=11.69, width=8.2)

