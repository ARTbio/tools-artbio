#!/usr/bin/env Rscript

# Description:
#   This script serves as the backend for the Galaxy FACETS Analysis tool.
#   It takes a SNP pileup file as input and performs allele-specific copy
#   number analysis using the R package 'facets'.
# ==============================================================================

# --- Load Libraries ---
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(facets))

# --- Source the external plot_facets_enhanced function ---
# This finds the path of the currently running script and sources
# the R function file relative to it.
initial_opts <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=", "", initial_opts[grep("--file=", initial_opts)]))
source(file.path(script_path, "plot_facets_enhanced-v22.R"))

# --- Define and Parse Arguments ---

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
parser$add_argument("--output_vcf",
    type = "character", required = TRUE,
    help = "Path for the output VCF file with CNV calls."
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
parser$add_argument("--enable_merging",
    action = "store_true", default = FALSE,
    help = "If specified, enables the post-processing step to merge adjacent and similar CNV segments."
)
parser$add_argument("--merge_gap_abs",
    type = "integer", default = 1000000,
    help = "Absolute maximum gap in bp to merge adjacent CNV segments."
)
parser$add_argument("--merge_gap_rel",
    type = "double", default = 0.5,
    help = "Relative maximum gap (fraction of avg. segment length) to merge segments."
)
parser$add_argument("--vcf_min_nhet",
    type = "integer", default = 2,
    help = "VCF Post-Filter: Minimum number of heterozygous SNPs for a segment to be kept."
)
parser$add_argument("--vcf_min_num_mark",
    type = "integer", default = 3,
    help = "VCF Post-Filter: Minimum number of total markers for a segment to be kept."
)
#' Classify CNV segments based on TCN/LCN
classify_cnv <- function(cncf_df) {
    cncf_df$sv_type <- NA_character_
    cncf_df$sv_type[cncf_df$tcn.em == 2 & (cncf_df$lcn.em == 1 | is.na(cncf_df$lcn.em))] <- "NEUTR"
    cncf_df$sv_type[is.na(cncf_df$sv_type) & cncf_df$tcn.em > 2] <- "DUP"
    cncf_df$sv_type[is.na(cncf_df$sv_type) & cncf_df$tcn.em < 2 & !is.na(cncf_df$lcn.em) & cncf_df$lcn.em > 0] <- "HEMIZYG_DEL"
    cncf_df$sv_type[is.na(cncf_df$sv_type) & cncf_df$tcn.em < 2 & !is.na(cncf_df$lcn.em) & cncf_df$lcn.em == 0] <- "HOMOZYG_DEL"
    cncf_df$sv_type[is.na(cncf_df$sv_type) & cncf_df$tcn.em == 2 & !is.na(cncf_df$lcn.em) & cncf_df$lcn.em == 0] <- "CN_LOH"

    # Remplacer les NA restants (si tcn.em < 2 mais lcn.em est NA) par un type général
    cncf_df$sv_type[is.na(cncf_df$sv_type) & cncf_df$tcn.em < 2] <- "DEL"

    return(cncf_df)
}

#' Create a VCF header (explicit version)
create_vcf_header <- function(sample_id, purity, ploidy) {
    header <- c(
        "##fileformat=VCFv4.2",
        paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
        paste0("##source=FACETS_v", packageVersion("facets")),
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant (standard VCF tags: DEL, DUP, CNV)\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">",
        # --- MODIFICATION ---
        "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"FACETS event classification. Possible values: DUP, HEMIZYG_DEL, HOMOZYG_DEL, CN_LOH\">",
        # --- FIN MODIFICATION ---
        "##INFO=<ID=TCN,Number=1,Type=Integer,Description=\"Total Copy Number (EM fit)\">",
        "##INFO=<ID=LCN,Number=1,Type=Integer,Description=\"Lesser Copy Number (EM fit)\">",
        "##INFO=<ID=NUM_MARK,Number=1,Type=Integer,Description=\"Number of SNPs in the segment\">",
        "##INFO=<ID=NHET,Number=1,Type=Integer,Description=\"Number of heterozygous SNPs in the segment\">",
        paste0("##FACETS_PURITY=", round(purity, 4)),
        paste0("##FACETS_PLOIDY=", round(ploidy, 4)),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    return(header)
}

# ==============================================================================
# Function to merge adjacent and similar CNV segments
#
# This function implements a merging algorithm that reflects human
# by using a hybrid proximity condition. Two segments are merged
# if they have the same CNV state and are close to each other, both in absolute
# and relative terms.
#
# @param cnv_df A data frame of CNV calls, expected to have columns like
#   'chrom', 'start', 'end', 'svtype', 'tcn.em', 'lcn.em', etc.
# @param max_gap_abs An integer. The absolute maximum distance (in bp) allowed
#   between two segments to consider them for merging. Acts as a safeguard.
# @param max_gap_rel A numeric value (0 to 1). The maximum relative distance,
#   expressed as a fraction of the average length of the two adjacent segments.
# @return A new data frame with the similar adjacent segments merged.
# ==============================================================================
merge_segments <- function(cnv_df, max_gap_abs = 1000000, max_gap_rel = 0.5) {
    # If there's nothing to merge, return the original data frame
    if (nrow(cnv_df) < 2) {
        return(cnv_df)
    }

    # Ensure the data frame is sorted by genomic position
    cnv_df <- cnv_df[order(cnv_df$chrom, cnv_df$start), ]

    merged_rows <- list()
    current_row <- cnv_df[1, ]

    for (i in 2:nrow(cnv_df)) {
        next_row <- cnv_df[i, ]

        # Basic criteria: segments must be of the same type and CN state
        same_chrom <- current_row$chrom == next_row$chrom
        same_svtype <- current_row$svtype == next_row$svtype
        same_event <- current_row$event == next_row$event
        same_tcn <- current_row$tcn.em == next_row$tcn.em
        same_lcn <- identical(current_row$lcn.em, next_row$lcn.em) # Handles NA safely

        # If basic criteria are met, evaluate proximity
        if (same_chrom && same_svtype && same_event && same_tcn && same_lcn) {
            gap <- next_row$start - current_row$end

            # Calculate the relative threshold based on the average size of the two segments
            len_a <- current_row$end - current_row$start
            len_b <- next_row$end - next_row$start
            relative_threshold <- ((len_a + len_b) / 2) * max_gap_rel

            # Hybrid merge condition: gap must be below BOTH absolute and relative thresholds
            if (gap >= 0 && gap <= max_gap_abs && gap <= relative_threshold) {
                # Merge: update the end of the current segment and aggregate numeric fields
                current_row$end <- next_row$end
                current_row$num.mark <- current_row$num.mark + next_row$num.mark
                current_row$nhet <- current_row$nhet + next_row$nhet

                # Skip to the next iteration to try merging the newly expanded segment
                # with the one that follows.
                next
            }
        }

        # If no merge occurred, the current segment is final. Save it.
        merged_rows <- append(merged_rows, list(current_row))
        # The next segment becomes the new current segment.
        current_row <- next_row
    }

    # Append the very last segment (which is either a standalone or the result of the last merge)
    merged_rows <- append(merged_rows, list(current_row))

    # Reconstruct a single data frame from the list of merged rows
    do.call(rbind, merged_rows)
}
# --- Main Analysis Function ---
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
    plot_facets_enhanced(oo, emfit = fit, plot.type = "em", sname = args$sample_id)
    dev.off()
    png(file = args$output_spider, width = 8, height = 8, units = "in", res = 300)
    logRlogORspider(oo$out, oo$dipLogR)
    dev.off()

    # --- Generate VCF file ---

    # Classify segments and define standard SVTYPEs + detailed EVENTs
    cncf_for_vcf <- fit$cncf
    if (nrow(cncf_for_vcf) > 0) {
        cncf_for_vcf$svtype <- NA_character_
        cncf_for_vcf$event <- NA_character_

        # Duplications
        cncf_for_vcf[cncf_for_vcf$tcn.em > 2, c("svtype", "event")] <- c("DUP", "DUP")

        # Deletions
        cncf_for_vcf[cncf_for_vcf$tcn.em < 2, c("svtype")] <- "DEL"
        cncf_for_vcf[cncf_for_vcf$tcn.em == 1, c("event")] <- "HEMIZYG_DEL"
        cncf_for_vcf[cncf_for_vcf$tcn.em == 0, c("event")] <- "HOMOZYG_DEL"

        # Copy-Neutral LOH
        cncf_for_vcf[cncf_for_vcf$tcn.em == 2 & !is.na(cncf_for_vcf$lcn.em) & cncf_for_vcf$lcn.em == 0, c("svtype", "event")] <- c("CNV", "CN_LOH")

        # Filter normal segments (where'svtype' is still NA)
        cnv_calls <- cncf_for_vcf[!is.na(cncf_for_vcf$svtype), ]
    } else {
        cnv_calls <- data.frame()
    }


    if (nrow(cnv_calls) > 0) {
        if (args$enable_merging) {
            cnv_calls <- merge_segments(
                cnv_df = cnv_calls,
                max_gap_abs = args$merge_gap_abs,
                max_gap_rel = args$merge_gap_rel
            )
        }
        # Apply VCF post-filters to remove low-quality/artefactual segments
        # This addresses the issue of FACETS' EM algorithm sometimes creating
        # micro-segments that bypass the initial min.nhet segmentation parameter.
        original_rows <- nrow(cnv_calls)
        cnv_calls <- cnv_calls[
            cnv_calls$nhet >= args$vcf_min_nhet &
                cnv_calls$num.mark >= args$vcf_min_num_mark,
        ]
        cat(paste("Applied VCF post-filters: kept", nrow(cnv_calls), "of", original_rows, "segments.\n"))

        vcf_header <- create_vcf_header(args$sample_id, fit$purity, fit$ploidy)

        vcf_body <- apply(cnv_calls, 1, function(seg) {
            cnv_calls <- merge_segments(cnv_calls)
            alt_allele <- paste0("<", seg["svtype"], ">")
            info <- paste0(
                "END=", seg["end"],
                ";SVTYPE=", seg["svtype"],
                ";SVLEN=", as.integer(seg["end"]) - as.integer(seg["start"]),
                ";TCN=", seg["tcn.em"],
                ";LCN=", ifelse(is.na(seg["lcn.em"]), ".", seg["lcn.em"]),
                ";EVENT=", seg["event"],
                ";NUM_MARK=", seg["num.mark"],
                ";NHET=", seg["nhet"]
            )
            # Remove any space(s) immediately following an '=' sign in the INFO string.
            info <- gsub("=\\s+", "=", info)

            paste(seg["chrom"], seg["start"], ".", "N", alt_allele, ".", "PASS", info, sep = "\t")
        })

        writeLines(c(vcf_header, vcf_body), con = args$output_vcf)
    } else {
        vcf_header <- create_vcf_header(args$sample_id, fit$purity, fit$ploidy)
        writeLines(vcf_header, con = args$output_vcf)
    }
}

# --- Execution Block ---
if (!interactive()) {
    args <- parser$parse_args()
    main(args)
}
