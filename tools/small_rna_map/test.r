# Table is the data frame
p <- ggplot(Table, aes(x=`Nbr_reads`, y=Coordinate, label=Nbr_reads)) + 
  geom_point(stat='identity', fill="black", size=6)  +
  geom_segment(aes(y = 0, 
                   x = `Nbr_reads`, 
                   yend = Coordinate, 
                   xend = `Nbr_reads`), 
               color = "black") +
  geom_text(color="white", size=2) +
  facet_wrap(~Chromosome, scales="free")+
  coord_flip()

p


