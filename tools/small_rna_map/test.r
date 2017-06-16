# Table is the data frame

library("ggplot2")
library("gridExtra")
library("RColorBrewer")
library("gtable")
library("grid")

theme_set(theme_bw())
#Table=read.delim(your_input, header=T, row.names=NULL)
Table <- within(Table[1:27,], Nbr_reads[Polarity=="R"] <- (Nbr_reads[Polarity=="R"]*-1)) 

myColors <- brewer.pal(3,"Set1")
names(myColors) <- levels(Table$Polarity)
colScale <- scale_colour_manual(name = "Polarity",values = myColors)

#Make initial figures
p1 <- ggplot(Table, aes(x=Coordinate, y=Nbr_reads, colour=Polarity)) +
  colScale+
  geom_segment(aes(y = 0, 
                   x = Coordinate, 
                   yend = Nbr_reads, 
                   xend = Coordinate,
			 color=Polarity),
			 alpha=1
               ) +
  geom_segment(aes(y = -Nbr_reads, 
			 x = 0,
			 yend=Nbr_reads,
			 xend=Chrom_length), alpha=0
		   )+
  facet_wrap(Dataset~Chromosome, scales="free")+
  geom_hline(yintercept=0, size=0.3)#+

p2 <- ggplot(Table, aes(x = Coordinate, y = Median)) +
  geom_point(aes(col=Median),size = 1, colour="black") + 
  expand_limits(y = seq(0,max(Table$Median),by=4)) +
  facet_wrap(Dataset~Chromosome, scales="free") +
  geom_segment(aes(y = -Nbr_reads, 
			 x = 0,
			 yend=Nbr_reads,
			 xend=Chrom_length), alpha=0
		   )+
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        legend.position = "bottom")

# Get the ggplot grobs
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# Get the locations of the plot panels in g1.
pp <- c(subset(g1$layout, grepl("panel", g1$layout$name), se = t:r))

# https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
# Overlap panels for second plot on those of the first plot
g <- gtable_add_grob(g1, g2$grobs[grepl("panel", g1$layout$name)], 
      pp$t, pp$l, pp$b, pp$l)


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
# Get the y axis title from g2
index <- which(g2$layout$name == "ylab-l") # Which grob contains the y axis title?   EDIT HERE
ylab <- g2$grobs[[index]]                # Extract that grob
ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications

# Put the transformed label on the right side of g1
g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, ylab, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, clip = "off", name = "ylab-r")

# Get the y axis from g2 (axis line, tick marks, and tick mark labels)
index <- which(g2$layout$name == "axis-l-1-1")  # Which grob.    EDIT HERE
yaxis <- g2$grobs[[index]]                    # Extract the grob

#swap tick marks and tick mark labels
ticks <- yaxis$children[[2]]
ticks$widths <- rev(ticks$widths)
ticks$grobs <- rev(ticks$grobs)

#move the tick marks

plot_theme <- function(p) {
  plyr::defaults(p$theme, theme_get())
}
tml <- plot_theme(p1)$axis.ticks.length   # Tick mark length

# Put ticks back into yaxis
yaxis$children[[2]] <- ticks

# Put the transformed yaxis on the right side of g1
g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, yaxis, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, 
    name = "axis-r")


grid.draw(g)


pdf(file="output.pdf", paper="special", height=11.69, width=8.2677)
p<-grid.draw(g)
dev.off()



