if (length(commandArgs(TRUE)) == 0) {
  system("Rscript barplot.r -h", intern = F)
  q("no")
}


# load packages that are provided in the conda env
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()
library(optparse)
library(ggplot2)
library(scales)



#Arguments
option_list = list(
  make_option(
    c("-i", "--input"),
    default = NA,
    type = 'character',
    help = "Input file that contains count data (no header)"
  ),
  make_option(
    c("-o", "--barplot"),
    default = NA,
    type = 'character',
    help = "PDF output file"
  )
)

opt = parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))


## 
annotations = read.delim(opt$input, header=F)
colnames(annotations) = c("class", "counts")
annotations = cbind(annotations, fraction=annotations$counts/annotations$counts[1])
annotations = annotations[-1,]
# ggplot2 plotting
ggplot(annotations, aes(x="classes", y=fraction, fill=class)) +
geom_bar(width = .7, position=position_stack(), stat = "identity") +
geom_text(aes(label = percent(fraction)), position = position_stack(vjust = 0.5),size = 4)
ggtitle('Class proportions') 
ggsave(file=opt$barplot, device="pdf")
