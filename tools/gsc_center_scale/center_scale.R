options(show.error.messages = FALSE,
        error = function(){
            cat(geterrmessage(), file=stderr())
            q("no", 1, FALSE)
            }
)
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)

# Arguments
option_list <- list(
  make_option(
    '--data',
    default = NA,
    type = 'character',
    help = "Input file that contains values to transform. Must be tabular separated,
            with columns and row names, variables in rows, observations in columns  [default : '%default' ]"
  ),
  make_option(
   '--center',
    default = TRUE,
    type = 'logical',
    help = "center data to the mean [default : '%default' ]"
  ),
  make_option(
   '--scale',
    default = TRUE,
    type = 'logical',
    help = "scale data to standard deviation [default : '%default' ]"
  ),
  make_option(
    '--factor',
    default = '',
    type = 'character',
    help = "A two-column observations|factor_levels table, to group observations to be transformed by levels  [default : '%default' ]"
  ),
  make_option(
    '--output',
    default = 'res.tab',
    type = 'character',
    help = "Table of transformed values [default : '%default' ]"
  )
)

transform <- function(df, center = TRUE, scale = TRUE) {
    transfo <- scale(
        t(df),
        center = center,
        scale = scale
        )
    return(as.data.frame(t(transfo)))
}

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

data <- read.delim(
    opt$data,
    check.names = FALSE,
    header = TRUE,
    row.names = 1,
    sep = '\t'
)

if (opt$factor != '') {
    data.factor = read.delim(
        opt$factor,
        check.names = FALSE,
        header = TRUE,
        sep = '\t'
        )
    colnames(data.factor) <- c("cellid", "level")
    data.transformed <- data.frame(row.names=rownames(data), stringsAsFactors = FALSE)
    for (group in levels(data.factor$level)){
        subcells <- as.data.frame(subset(data.factor, level == group, select = cellid))
        subdata <- as.data.frame(subset(data, select = as.vector(subcells$cellid)))
        subdata.transformed <- transform(subdata, center = as.logical(opt$center),
                                                  scale = as.logical(opt$scale))
        data.transformed <- cbind(data.transformed, subdata.transformed)
    }
} else {
    data.transformed <- transform(data, center = as.logical(opt$center),
                                        scale = as.logical(opt$scale))
}


write.table(
  cbind(gene = rownames(data.transformed), data.transformed),
  opt$output,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
