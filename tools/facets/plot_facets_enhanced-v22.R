plot_facets_enhanced <- function(x, emfit = NULL, clustered = FALSE, 
                                            plot.type = c("em", "naive", "both", "none"), 
                                            sname = NULL, add.legend = TRUE) {
  
  def.par <- par(no.readonly = TRUE)
  plot.type <- match.arg(plot.type)
  
  # Setup layout
  if (plot.type == "none") 
    layout(matrix(1:2, ncol = 1))
  if (plot.type == "em") 
    layout(matrix(rep(1:4, c(9, 9, 6, 1)), ncol = 1))
  if (plot.type == "naive") 
    layout(matrix(rep(1:4, c(9, 9, 6, 1)), ncol = 1))
  if (plot.type == "both") 
    layout(matrix(rep(1:6, c(9, 9, 6, 1, 6, 1)), ncol = 1))
  
  par(mar = c(0.25, 3, 0.25, 1), mgp = c(1.75, 0.6, 0), 
      oma = c(3, 0, 1.25, 0))
  
  # Prepare data
  jseg <- x$jointseg
  chrbdry <- which(diff(jseg$chrom) != 0)
  
  if (missing(emfit)) {
    out <- x$out
    if (plot.type == "em" | plot.type == "both") {
      warning("emfit is missing; plot.type set to naive")
      plot.type <- "naive"
    }
  } else {
    out <- emfit$cncf
    out$tcn <- x$out$tcn
    out$lcn <- x$out$lcn
    out$cf <- x$out$cf
  }
  
  if (clustered) {
    cnlr.median <- out$cnlr.median.clust
    mafR <- out$mafR.clust
    mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
  } else {
    cnlr.median <- out$cnlr.median
    mafR <- out$mafR
  }
  
  mafR <- abs(mafR)
  chrcol <- 1 + rep(out$chrom - 2 * floor(out$chrom/2), out$num.mark)
  nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
  segbdry <- cumsum(c(0, out$num.mark))
  segstart <- segbdry[-length(segbdry)]
  segend <- segbdry[-1]
  
  # Plot 1: log-ratio
  plot(jseg$cnlr[is.finite(jseg$cnlr)], pch = ".", cex = 2, 
       col = c("grey", "lightblue", "azure4", "slateblue")[chrcol], 
       ylab = "log-ratio", xaxt = "n")
  abline(v = chrbdry, lwd = 0.25)
  
  # Calculate reference lines
  global_median <- median(jseg$cnlr, na.rm = TRUE)
  diploid_logr <- x$dipLogR
  lines_diff <- abs(global_median - diploid_logr)
  
  # Only show reference lines if they are sufficiently different (> 0.1)
  legend_items <- c("Segments")
  legend_cols <- c(2)
  legend_lwd <- c(1.75)
  legend_lty <- c(1)
  
  if (lines_diff > 0.1) {
    # Lines are different enough to show both
    abline(h = global_median, col = "green2", lwd = 2.5)
    abline(h = diploid_logr, col = "magenta4", lwd = 2.5, lty = 2)
    legend_items <- c(legend_items, 
                      sprintf("Global median (%.3f)", global_median),
                      sprintf("Diploid LogR (%.3f)", diploid_logr))
    legend_cols <- c(legend_cols, "green2", "magenta4")
    legend_lwd <- c(legend_lwd, 2.5, 2.5)
    legend_lty <- c(legend_lty, 1, 2)
  } else {
    # Lines are too close, show only one with combined label
    abline(h = global_median, col = "green2", lwd = 2.5)
    legend_items <- c(legend_items, 
                      sprintf("Median/Diploid (%.3f)", global_median))
    legend_cols <- c(legend_cols, "green2")
    legend_lwd <- c(legend_lwd, 2.5)
    legend_lty <- c(legend_lty, 1)
  }
  
  segments(segstart, cnlr.median, segend, cnlr.median, lwd = 1.75, col = 2)
  
  # Add legend for log-ratio plot
  if (add.legend) {
    legend("topright", 
           legend = legend_items,
           col = legend_cols, 
           lwd = legend_lwd,
           lty = legend_lty,
           bty = "n", 
           cex = 0.7)
  }
  
  # Plot 2: log-odds-ratio
  plot(jseg$valor[is.finite(jseg$cnlr)], pch = ".", cex = 2.5, 
       col = c("grey", "lightblue", "azure4", "slateblue")[chrcol], 
       ylab = "log-odds-ratio", ylim = c(-4, 4), xaxt = "n")
  abline(v = chrbdry, lwd = 0.25)
  segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd = 1.75, col = 2)
  segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd = 1.75, col = 2)
  
  # Add legend for log-odds-ratio plot
  if (add.legend) {
    legend("topright", 
           legend = c("BAF segments"), 
           col = c(2), 
           lwd = c(1.75), 
           bty = "n", 
           cex = 0.7)
  }
  
  cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10), "bisque2")
  
  # Plot 3: copy number (naive)
  if (plot.type == "naive" | plot.type == "both") {
    out$tcn[out$tcn > 10] <- 9 + log10(out$tcn[out$tcn > 10])
    ii <- which(out$lcn > 5)
    if (length(ii) > 0) 
      out$lcn[ii] <- 5 + log10(out$lcn[ii])
    
    plot(c(0, length(jseg$cnlr)), c(0, max(out$tcn)), type = "n", 
         ylab = "copy number (nv)", xaxt = "n")
    abline(v = chrbdry, lwd = 0.25)
    segments(segstart, out$lcn, segend, out$lcn, lwd = 1.75, col = 2)
    segments(segstart, out$tcn, segend, out$tcn, lwd = 1.75, col = 1)
    
    # Add legends for copy number plot (including CF legend)
    if (add.legend) {
      # CN legend on top left
      legend("topleft", 
             legend = c("Total CN (TCN)", "Minor CN (LCN)"), 
             col = c(1, 2), 
             lwd = 1.75, 
             bty = "n", 
             cex = 0.7)
      
      # CF legend on top right
      legend("topright", 
             legend = c("CF = 0 (normal)", "CF = 0.2 (low tumor)", 
                        "CF = 0.8 (high tumor)", "CF = 1.0 (pure tumor)"), 
             fill = c("white", "lightblue", "steelblue", "tan"),
             border = "black",
             bty = "o",
             box.col = "black",
             box.lwd = 1.5,
             cex = 0.65,
             title = "Cellular Fraction (CF)",
             bg = "white")
    }
    
    # CF bar (naive)
    plot(c(0, length(jseg$cnlr)), 0:1, type = "n", ylab = "", 
         xaxt = "n", yaxt = "n")
    mtext("cf-nv", side = 2, at = 0.5, line = 0.3, las = 2, cex = 0.75)
    cfcol <- cfpalette[round(10 * out$cf + 0.501)]
    rect(segstart, 0, segend, 1, col = cfcol, border = NA)
  }
  
  # Plot 4: copy number (EM)
  if (plot.type == "em" | plot.type == "both") {
    out$tcn.em[out$tcn.em > 10] <- 9 + log10(out$tcn.em[out$tcn.em > 10])
    ii <- which(out$lcn.em > 5)
    if (length(ii) > 0) 
      out$lcn.em[ii] <- 5 + log10(out$lcn.em[ii])
    
    plot(c(0, length(jseg$cnlr)), c(0, max(out$tcn.em)), 
         type = "n", ylab = "copy number (em)", xaxt = "n")
    abline(v = chrbdry, lwd = 0.25)
    segments(segstart, out$lcn.em, segend, out$lcn.em, lwd = 1.75, col = 2)
    segments(segstart, out$tcn.em, segend, out$tcn.em, lwd = 1.75, col = 1)
    
    # Add legends for EM copy number plot (including CF legend)
    if (add.legend) {
      # CN legend on top left
      legend("topleft", 
             legend = c("Total CN (TCN)", "Minor CN (LCN)"), 
             col = c(1, 2), 
             lwd = 1.75, 
             bty = "n", 
             cex = 0.7)
      
      # CF legend on top right
      legend("topright", 
             legend = c("CF = 0 (normal)", "CF = 0.2 (low tumor)", 
                        "CF = 0.8 (high tumor)", "CF = 1.0 (pure tumor)"), 
             fill = c("white", "lightblue", "steelblue", "tan"),
             border = "black",
             bty = "o",
             box.col = "black",
             box.lwd = 1.5,
             cex = 0.65,
             title = "Cellular Fraction (CF)",
             bg = "white")
    }
    
    # CF bar (EM)
    plot(c(0, length(jseg$cnlr)), 0:1, type = "n", ylab = "", 
         xaxt = "n", yaxt = "n")
    mtext("cf-em", side = 2, at = 0.5, line = 0.2, las = 2, cex = 0.75)
    cfcol <- cfpalette[round(10 * out$cf.em + 0.501)]
    rect(segstart, 0, segend, 1, col = cfcol, border = NA)
  }
  
  # X-axis with chromosome labels
  chromlevels <- x$chromlevels
  if (is.null(chromlevels)) 
    chromlevels <- 1:length(nn)
  axis(labels = chromlevels, side = 1, 
       at = (nn + c(0, nn[-length(nn)]))/2, cex = 0.65)
  mtext(side = 1, line = 1.75, "Chromosome", cex = 0.8)
  
  # Sample name title
  if (!missing(sname)) 
    mtext(sname, side = 3, line = 0, outer = TRUE, cex = 0.8)
  
  par(def.par)
}
