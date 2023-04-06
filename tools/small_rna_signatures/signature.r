## Setup R error handling to go to stderr
options(show.error.messages = FALSE,
        error = function() {
            cat(geterrmessage(), file = stderr())
            q("no", 1, FALSE)
        }
)
warnings()

library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(optparse)

option_list <- list(
    make_option("--h_dataframe", type = "character",
                help = "path to h-signature dataframe"),
    make_option("--z_dataframe", type = "character",
                help = "path to z-signature dataframe"),
    make_option("--plot_method", type = "character",
                help = "How  data should be plotted (global or lattice)"),
    make_option("--pdf", type = "character", help = "path to the pdf file with plots"),
    make_option("--title", type = "character", help = "Graph Title"),

    make_option("--npairs_ylim", type = "integer", default = 0, help = "Maximum of Y-scale of the Numbers-of-pairs plot"),
    make_option("--npairszscore_ylim", type = "integer", default = 0, help = "Maximum of Y-scale of the Number-of-pairs-Z−scores plot"),
    make_option("--overlapprob_ylim", type = "integer", default = 0, help = "Maximum of Y-scale of the Overlap-probabilities plot"),
    make_option("--overlapprobzscore_ylim", type = "integer", default = 0, help = "Maximum of Y-scale of the Overlap-probabilities-Z−scores plot")
    )
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# data frames implementation
h_dataframe <- read.delim(args$h_dataframe, header = FALSE)
colnames(h_dataframe) <- c("chrom", "overlap", "sig", "z-score")
h_dataframe$sig <- h_dataframe$sig * 100  # to get probs in %
z_dataframe <- read.delim(args$z_dataframe, header = FALSE)
colnames(z_dataframe) <- c("chrom", "overlap", "sig", "z-score")

# graph limits parameters
if (args$npairs_ylim == 0) {
    npairs_ylim <- NULL
    } else {
    npairs_ylim <- c(0, args$npairs_ylim)
}

if (args$npairszscore_ylim == 0) {
    npairszscore_ylim <- NULL
    } else {
    npairszscore_ylim <- c(0, args$npairszscore_ylim)
}

if (args$overlapprob_ylim == 0) {
    overlapprob_ylim <- NULL
    } else {
    overlapprob_ylim <- c(0, args$overlapprob_ylim)
}

if (args$overlapprobzscore_ylim == 0) {
    overlapprobzscore_ylim <- NULL
    } else {
    overlapprobzscore_ylim <- c(0, args$overlapprobzscore_ylim)
}


# functions
    globalgraph <- function() {
        pdf(args$pdf)
        par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))

            plot(z_dataframe[z_dataframe$chrom == "all_chromosomes", c(2, 3)],
                     type = "h", main = "Numbers of pairs", cex.main = 1, xlab = "overlap (nt)",
                     ylab = "Numbers of pairs", col = "darkslateblue", lwd = 4,
                     ylim = npairs_ylim
                     )

            plot(z_dataframe[z_dataframe$chrom == "all_chromosomes", c(2, 4)],
                     type = "l", main = "Number of pairs Z-scores", cex.main = 1, xlab = "overlap (nt)",
                     ylab = "z-score", pch = 19, cex = 0.2, col = "darkslateblue", lwd = 2,
                     ylim = npairszscore_ylim
                     )

            plot(h_dataframe[h_dataframe$chrom == "all_chromosomes", c(2, 3)],
                     type = "l", main = "Overlap probabilities", cex.main = 1,
                     xlab = "overlap (nt)",
                     ylab = "Probability [%]", pch = 19,
                     col = "darkslateblue", lwd = 2,
                     ylim = overlapprob_ylim
                     )

            plot(h_dataframe[h_dataframe$chrom == "all_chromosomes", c(2, 4)],
                     type = "l", main = "Overlap Probability Z-scores", cex.main = 1,
                     xlab = "overlap (nt)", ylab = "z-score", pch = 19, cex = 0.2,
                     col = "darkslateblue", lwd = 2,
                     ylim = overlapprobzscore_ylim
                     )

        mtext(args$title, outer = TRUE, cex = 1)
        dev.off()
    }

    treillisgraph <- function(df, ...) {
          pdf(args$pdf, paper = "special", height = 11.69, width = 6)
          p <- xyplot(sig ~ overlap | factor(method, levels = unique(method)) + chrom,
                   data = df,
                   type = "l",
                   col = "darkblue",
                   cex = 0.5,
                   scales = list(y = list(tick.number = 4, relation = "free", cex = 0.6,
                                          rot = 0),
                                 x = list(cex = 0.6, alternating = FALSE)),
                   xlab = "Overlap",
                   ylab = "signature (Nbr of pairs / Overlap prob.)",
                   main = args$title,
                   par.strip.text = list(cex = .5),
                   pch = 19, lwd = 2,
                   as.table = TRUE,
                   layout = c(2, 12),
                   newpage = TRUE,
                   ...)
           plot(p)
           dev.off()
    }

# main

if (args$plot_method == "global") {
    globalgraph()
}

if (args$plot_method == "lattice") {
    # rearrange dataframes
    h_sig <- h_dataframe[, c(1, 2, 3)]
    h_sig <- cbind(rep("Overlap Prob (%)", length(h_sig[, 1])), h_sig)
    colnames(h_sig) <- c("method", "chrom", "overlap", "sig")
    z_pairs <- z_dataframe[, c(1, 2, 3)]
    z_pairs <- cbind(rep("Nbr of pairs", length(z_pairs[, 1])), z_pairs)
    colnames(z_pairs) <- c("method", "chrom", "overlap", "sig")
    lattice_df <- rbind(z_pairs, h_sig)
    par_settings_treillis <- list(strip.background = list(
            col = c("lightblue", "lightgreen")))
    treillisgraph(lattice_df, par.settings = par_settings_treillis)
}
