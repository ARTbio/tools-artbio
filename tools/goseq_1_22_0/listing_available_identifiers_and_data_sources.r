library("reshape2")
# install all packages
packages=read.table("org_packages.tab")
biocLite(packages$V1, dependencies=TRUE)
# load all packages
lapply(packages$V1, require, character.only = TRUE)
# create package_name vector
pkg_name = sapply(packages$V1, function(x) paste("package:", x, sep="") )
# list package_functions
list_package_f = sapply(pkg_name, ls)
organism_names = sapply(packages$V1, function(x) paste(eval( parse( text=paste( gsub (".db$", "", x ), "ORGANISM", sep="") ) ), paste( paste("(", x, sep=""), ")", sep="")))
identifiers = c( "GENENAME", "UNIGENE", "REFSEQ", "SYMBOL", "ENSEMBL" )
org_name_tab = data.frame(organism_names, packages$V1)
categories = c( "PMID","ENZYME", "GO" )

# get dataframe suitable for galaxy's <filter></> tagset
filter_tab = melt(list_package_f)
filter_tab$L1 = sapply( filter_tab$L1, function(x) gsub( "package:", "", x) )
#filter_tab$L1 = sapply( filter_tab$L1, function(x) gsub( ".db$", "", x ) )
patterns=paste(unique(sapply(filter_tab$L1, function(x) gsub( ".db$", "", x )) ), collapse="|")
filter_tab[,1] = gsub( patterns, "", filter_tab[,1] )

available_identifiers = subset(filter_tab, value %in% identifiers)
available_categories = subset(filter_tab, value %in% categories)

write.table(available_identifiers, file = "available_identifiers.tab", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(available_categories, file = "available_categories.tab", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(org_name_tab, file = "org_name.tab", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)