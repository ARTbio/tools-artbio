# Use this script to generate the .loc.sample with available organisms, available gene identifier and available categories
library("reshape2")
# install all packages
packages=read.table("org_packages.tab")
biocLite(packages$V1, dependencies=TRUE)
# load all packages
lapply(packages$V1, require, character.only = TRUE)
# create package_name vector
pkg_name = sapply(packages$V1, function(x) paste("package:", x, sep="") )
# list package_functions
organism_names = sapply(packages$V1, function(x) paste(eval( parse( text=paste( gsub (".db$", "", x ), "ORGANISM", sep="") ) ), paste( paste("(", x, sep=""), ")", sep="")))
identifiers = c( "GENENAME", "UNIGENE", "UNIPROT", "REFSEQ", "SYMBOL", "ENSEMBL", "FLYBASECG", "ACCNUM" )
org_name_tab = data.frame(organism_names, packages$V1)
categories = c( "PMID", "ENZYME", "GO2ALLEGS", "PATH", "GO2ALLTAIRS", "GO2ALLORFS", "PFAM", "PROSITE" )

# get dataframe suitable for galaxy's <filter></> tagset
filter_tab = melt(sapply(pkg_name, ls))
filter_tab$L1 = sapply( filter_tab$L1, function(x) gsub( "package:", "", x) )
patterns=paste(unique(sapply(filter_tab$L1, function(x) gsub( ".db$", "", x )) ), collapse="|")
filter_tab[,1] = gsub( patterns, "", filter_tab[,1] )

# add the ENTREZ id format to the available_identifiers

available_identifiers = subset(filter_tab, value %in% identifiers)
available_identifiers = cbind(available_identifiers[,1], available_identifiers)
available_categories = subset(filter_tab, value %in% categories)
available_categories = cbind(available_categories[,1], available_categories)
entrez = data.frame(rep("ENTREZ", length(packages$V1)), rep("ENTREZ", length(packages$V1)), packages$V1)
colnames(entrez) = colnames(available_identifiers)
available_identifiers = rbind(entrez, available_identifiers)

write.table(available_identifiers, file = "available_identifiers.loc.sample", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(available_categories, file = "available_categories.loc.sample", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(org_name_tab, file = "org_name.loc.sample", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)