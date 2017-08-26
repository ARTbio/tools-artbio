## Setup R error handling to go to stderr
options( show.error.messages=F,
       error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
warnings()

library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(optparse)

option_list <- list(
    make_option("--h_dataframe", type="character", help="path to h-signature dataframe"),
    make_option("--z_dataframe", type="character", help="path to z-signature dataframe"),
    make_option("--plot_method", type = "character", help="How  data should be plotted (global of by-item)"),
    make_option("--pdf", type = "character", help="path to the pdf file with plots"),
    make_option("--title", type = "character", help="Graph Title")
    )
 
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args = parse_args(parser)



 
# data frames implementation
h_dataframe = read.delim(args$h_dataframe, header=F)
colnames(h_dataframe) = c("chrom", "overlap", "prob", "z-score")
h_dataframe$prob = h_dataframe$prob * 100  # to get probs in %
z_dataframe = read.delim(args$z_dataframe, header=F)
colnames(z_dataframe) = c("chrom", "overlap", "nbre_pairs", "z-score")

# functions
      globalgraph = function () {
        pdf(args$pdf)
        par(mfrow=c(2,2),oma = c(0, 0, 3, 0))
        
        plot(z_dataframe[z_dataframe$chrom == "all_chromosomes", c(2,3)],
             type = "h", main="Numbers of pairs", cex.main=1, xlab="overlap (nt)",
             ylab="Numbers of pairs", col="darkslateblue", lwd=4)

        plot(z_dataframe[z_dataframe$chrom == "all_chromosomes", c(2,4)],
             type = "l", main="Number of pairs Z-scores", cex.main=1, xlab="overlap (nt)",
             ylab="z-score", pch=19, cex=0.2, col="darkslateblue", lwd=2)

        plot(h_dataframe[h_dataframe$chrom == "all_chromosomes", c(2,3)],
             type = "l", main="Overlap probabilities", cex.main=1, xlab="overlap (nt)",
             ylab="Probability [%]", ylim=c(0,50), pch=19, col="darkslateblue", lwd=2)

        plot(h_dataframe[h_dataframe$chrom == "all_chromosomes", c(2,4)],
             type = "l", main="Overlap Probability Z-scores", cex.main=1,
             xlab="overlap (nt)", ylab="z-score", pch=19, cex=0.2,
             col="darkslateblue", lwd=2)

        mtext(args$title, outer = TRUE, cex=1)
        dev.off()
      }

      treillisgraph = function () {
        pdf( args$pdf, paper="special", height=11.69, width=8.2677 )
#        signature = read.delim("${output}", header=TRUE)
#        xyplot(signature[,3]*100~signature[,1]|signature[,4], type = "l", xlim=c(${minscope},${maxscope}), main="ping-pong Signature of ${minquery}-${maxquery} against ${mintarget}-${maxtarget}nt small RNAs",
#             par.strip.text=list(cex=.5), strip=strip.custom(which.given=1, bg="lightblue"), scales=list(cex=0.5),
#             cex.main=1, cex=.5, xlab="overlap (nt)", ylab="ping-pong signal [%]",
#             pch=19, col="darkslateblue", lwd =1.5, cex.lab=1.2, cex.axis=1.2,
#             layout=c(4,12), as.table=TRUE, newpage = T)
#        dev.off()
#      }

      if (args$plot_method=="global") {
        globalgraph()

      }
      if(args$plot_method=="by-item") {
        treillisgraph()
      }

