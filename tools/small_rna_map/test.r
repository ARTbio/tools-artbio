# Table is the data frame
library("ggplot2")
library("gridExtra")
library(ggplus)

theme_set(theme_bw())
#Table=read.delim(your_input, header=T, row.names=NULL)
Table <- within(Table[1:27,], Nbr_reads[Polarity=="R"] <- (Nbr_reads[Polarity=="R"]*-1)) 
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
  geom_hline(yintercept=0, size=0.3)

p2 <- ggplot(Table, aes(x = Coordinate, y = Median)) +
  geom_point(aes(col=Median),size = 2, colour="black") + 
  expand_limits(y = seq(0,max(Table$Median),by=4)) +
  facet_wrap(Dataset~Chromosome, scales="free") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        legend.position = "bottom")


# Get the ggplot grobs
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

# Get the locations of the plot panels in g1.
pp <- c(subset(g1$layout, grepl("panel", g1$layout$name), se = t:r))

# Overlap panels for second plot on those of the first plot
g <- gtable_add_grob(g1, g2$grobs[grepl("panel", g1$layout$name)], 
      pp$t, pp$l, pp$b, pp$l)

pdf(file="output.pdf", paper="special", height=11.69, width=8.2677)
p<-grid.draw(g)
dev.off()



