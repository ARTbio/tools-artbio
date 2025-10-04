#!/usr/bin/env Rscript

# Description:
#   This script serves as the backend for the Galaxy FACETS Analysis tool.
#   It takes a SNP pileup file as input and performs allele-specific copy
#   number analysis using the R package 'facets'.
# ==============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(facets))

# --- 2. Define and Parse Arguments ---

# Create the parser
parser <- ArgumentParser(description = "Run FACETS algorithm on a SNP pileup file.")

# Define arguments
parser$add_argument("--pileup",
    type = "character", required = TRUE,
    help = "Path to the gzipped pileup CSV file."
)
parser$add_argument("--sample_id",
    type = "character", required = TRUE,
    help = "Sample ID used for plot titles and metadata."
)

parser$add_argument("--output_seg",
    type = "character", required = TRUE,
    help = "Path for the output segmentation file (TSV)."
)
parser$add_argument("--output_summary",
    type = "character", required = TRUE,
    help = "Path for the output summary file (TSV)."
)
parser$add_argument("--output_plots",
    type = "character", required = TRUE,
    help = "Path for the main output plots file (PNG)."
)
parser$add_argument("--output_spider",
    type = "character", required = TRUE,
    help = "Path for the diagnostic spider plot file (PNG)."
)
parser$add_argument("--cval",
    type = "double", default = 150,
    help = "Critical value for segmentation."
)
parser$add_argument("--min_nhet",
    type = "integer", default = 25,
    help = "Minimum number of heterozygous SNPs per segment."
)
parser$add_argument("--snp_nbhd",
    type = "integer", default = 300,
    help = "SNP neighborhood size for pre-processing. Crucial for sparse VCFs."
)
parser$add_argument("--gbuild",
    type = "character", default = "hg38",
    choices = c("hg19", "hg38", "hg18", "mm9", "mm10"),
    help = "Genome build used for alignment."
)

# --- 3. Main Analysis Function ---
main <- function(args) {
    # Set seed for reproducibility
    set.seed(1965)

    # --- Read the data with readSnpMatrix() from facets ---
    rcmat <- readSnpMatrix(args$pileup)

    # --- Pre-process sample ---
    xx <- preProcSample(rcmat, gbuild = args$gbuild, snp.nbhd = args$snp_nbhd)

    # --- Process sample (segmentation) ---
    oo <- procSample(xx, cval = args$cval, min.nhet = args$min_nhet)

    # --- Estimate ploidy/purity ---
    fit <- emcncf(oo)

    # Write the main segmentation file
    cncf_output <- fit$cncf
    if (nrow(cncf_output) > 0) {
        cncf_output$purity <- fit$purity
        cncf_output$ploidy <- fit$ploidy
        # Reorder columns to have purity/ploidy first for clarity
        cncf_output <- cncf_output[, c("purity", "ploidy", setdiff(names(cncf_output), c("purity", "ploidy")))]
    }
    write.table(cncf_output, file = args$output_seg, sep = "\t", quote = FALSE, row.names = FALSE)

    # Write a key-value summary file
    # A NULL value is replaced by NA to preserve vector length.
    summary_df <- data.frame(
        Parameter = c(
            "sample_id", "purity", "ploidy", "dipLogR", "loglik",
            "cval_param", "min_nhet_param", "snp_nbhd_param", "gbuild_param"
        ),
        Value = c(
            args$sample_id,
            ifelse(is.null(fit$purity), NA, fit$purity),
            ifelse(is.null(fit$ploidy), NA, fit$ploidy),
            ifelse(is.null(fit$dipLogR), NA, fit$dipLogR),
            ifelse(is.null(fit$loglik), NA, fit$loglik),
            args$cval,
            args$min_nhet,
            args$snp_nbhd,
            args$gbuild
        )
    )
    write.table(summary_df, file = args$output_summary, sep = "\t", quote = FALSE, row.names = FALSE)

    # Generate the plots PNG
    png(file = args$output_plots, width = 12, height = 8, units = "in", res = 300)
    plotSample(x = oo, emfit = fit, sname = args$sample_id)
    dev.off()
    png(file = args$output_spider, width = 8, height = 8, units = "in", res = 300)
    logRlogORspider(oo$out, oo$dipLogR)
    dev.off()

}

# --- 4. Execution Block ---
if (!interactive()) {
    args <- parser$parse_args()
    main(args)
}
