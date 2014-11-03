bioclite = "http://bioconductor.org/biocLite.R"
install.packages(c("lattice","latticeExtra"),dependencies=T,repos=http://cran.us.r-project.org)
source(bioclite)
installme=c(lattice,latticeExtra)
biocLite()
biocLite(installme)
quit(save="no")

