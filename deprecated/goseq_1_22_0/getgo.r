options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("goseq")
    library("optparse")
    library("reshape2")
})

sink(stdout(), type = "message")

option_list <- list(
    make_option(c("-p", "--package"), type="character", help = "Genome [used for looking up GO categories]"),
    make_option(c("-i", "--gene_id"), type="character", help="Gene ID format"),
    make_option(c("-c", "--cats"), type="character", help="Comma-seperated list of categories to fetch"),
    make_option(c("-o", "--output"), type="character", help="Path to output file")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# vars

package = args$package
gene_id = args$gene_id
output = args$output
cats = unlist(strsplit(args$cats, ','))

get_categories = function(package_str, gen, cat) {
  # gen should be ENSEMBL, UNIGENE, REFSEQ, SYMBOL or GENENAME
  # package should be org.Xx.eg.db
  # cat should be PMID, GO2ALLEGS, ENZYME or PATH
  library(package_str, character.only = TRUE)
  package = eval( parse( text=package_str ) )
  if( cat %in% c("GO2ALLEGS", "GO2ALLTAIRS", "GO2ALLORFS") ) {
    cat = "GOALL"
  }
  if(package_str == "org.Pf.plasmo.db") {
    keytype = "ORF"
    } else if(package_str == "org.At.tair.db") {
    keytype = "TAIR"
    } else {
    keytype = "ENTREZID"
    }
  entrez_cat = select(package, keys(package), cat, keytype)
  entrez_cat = entrez_cat[complete.cases(entrez_cat),]
    if( cat != "GOALL" ) {
      # add the origin of the term, so that there are no duplicate values e.g between ENZYME and PATH
      entrez_cat[,2] = sapply(entrez_cat[,2], function(x) paste(cat, x, sep=":"))
    } else {
      entrez_cat = entrez_cat[,c(1,2)] # we are discarding ontology (MF, CC, BP) and evidence class here
    }
  colnames(entrez_cat) = c(gen, "category")
  if( gen == "ENTREZ" ) {
    return( entrez_cat )
    } else {
      # We map ENTREZ to `gen`, but are potentially loosing gene identifiers where multiple identifiers match a single ENTREZ gene id.
      entrez_cat[,1] = mapIds(package, keys=as.character(entrez_cat[,1]), keytype=keytype, column=gen, multiVals="first")
      entrez_cat = entrez_cat[complete.cases(entrez_cat),]
      return(entrez_cat)
    }
}

result = lapply( cats, function(x) get_categories(package, gene_id, x ) )
result = do.call(rbind, result)

write.table(result, output, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
sessionInfo()
