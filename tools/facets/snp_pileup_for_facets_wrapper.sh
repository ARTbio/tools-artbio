#!/bin/bash

# ==============================================================================
#
# snp_pileup_for_facets_wrapper.sh (v4.0 - Optimal Parallelization)
#
# Description:
#   A highly optimized wrapper for snp-pileup. It divides the genome into
#   numerous, equally-sized chunks to ensure a balanced workload across all
#   allocated threads, maximizing performance and resource utilization.
#
# Author:
#   drosofff (ARTbio) & Gemini (Google)
#
# ==============================================================================

# --- Robustness settings ---
set -euo pipefail

# --- Help Function ---
usage() {
    cat << EOF
Usage: $0 -n <normal.bam> -t <tumor.bam> -v <snps.vcf.gz> -o <output.csv.gz> [OPTIONS]

REQUIRED ARGUMENTS:
  -n    Path to the normal BAM file (must be indexed).
  -t    Path to the tumor BAM file (must be indexed).
  -v    Path to the SNP reference VCF file (must be bgzip compressed and tabix indexed).
  -o    Path for the final output pileup file (will be bgzip compressed).

OPTIONS:
  -N    Number of parallel processes to use (default from \${GALAXY_SLOTS:-4}).
  -q    Minimum mapping quality for reads (default: 15).
  -Q    Minimum base quality for bases (default: 20).
  -p    Pseudo-SNP spacing in bp (default: 300).
  -A    Include anomalous read pairs (flag, default: not set).
  -h    Display this help message and exit.
EOF
    exit 0
}

# --- Argument Parsing ---
nprocs=${GALAXY_SLOTS:-4} # Default to Galaxy's variable, with a fallback
mapq=15
baseq=20
pseudo_snps=300
count_orphans=""

while getopts ":hn:t:v:o:N:q:Q:p:A" opt; do
    case ${opt} in
        h ) usage ;;
        n ) normal_bam="$OPTARG" ;;
        t ) tumor_bam="$OPTARG" ;;
        v ) snp_vcf="$OPTARG" ;;
        o ) output_pileup="$OPTARG" ;;
        N ) nprocs="$OPTARG" ;;
        q ) mapq="$OPTARG" ;;
        Q ) baseq="$OPTARG" ;;
        p ) pseudo_snps="$OPTARG" ;;
        A ) count_orphans="-A" ;;
        \? ) echo "Invalid option: -$OPTARG" >&2; usage ;;
        : ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# --- Tool path discovery ---
# (Assumes tools are in PATH)
SNP_PILEUP_EXE="snp-pileup"
SAMTOOLS_EXE="samtools"
BCFTOOLS_EXE="bcftools"

# --- Temp directory setup and cleanup trap ---
TMPDIR=$(mktemp -d)
trap 'rm -rf -- "$TMPDIR"' EXIT

# ==============================================================================
# --- Phase 1: Parallelization Strategy using Balanced Regions ---
# ==============================================================================

echo "Generating balanced genomic chunks for parallel processing..."

# The number of chunks is a multiple of the number of threads for optimal balance.
NUM_CHUNKS=$((nprocs * 10))

# Step 1: Get the list of chromosomes to process directly from the VCF data body.
# This is the most robust method. We then apply a positive filter to keep only
# autosomes and sex chromosomes.
VCF_PRIMARY_CHROMS=$( \
    zcat "${snp_vcf}" \
    | grep -v "^#" \
    | cut -f 1 \
    | sort -u \
    | grep -E '^(chr)?([0-9]+|X|Y)$' \
    || true \
)

if [ -z "$VCF_PRIMARY_CHROMS" ]; then
    echo "Error: No primary autosomes or sex chromosomes (1-22, X, Y) were found in the VCF file body." >&2
    exit 1
fi
echo "Found the following chromosomes to process:"
echo "$VCF_PRIMARY_CHROMS"

# Step 2: Get the geometry (lengths) for these specific chromosomes from the BAM header.
GENOME_GEOMETRY=$( \
    "${SAMTOOLS_EXE}" view -H "$normal_bam" \
    | awk -F'\t' '/^@SQ/ {print $2"\t"$3}' \
    | sed 's/SN://' | sed 's/LN://' \
    | grep -wFf <(echo "$VCF_PRIMARY_CHROMS") \
)

# Step 3: Calculate the total genome size to process and the "ideal" chunk size.
TOTAL_SIZE=$(echo "$GENOME_GEOMETRY" | awk '{sum+=$2} END {print sum}')
if [[ -z "$TOTAL_SIZE" || "$TOTAL_SIZE" -eq 0 ]]; then
    echo "Error: Could not determine genome size. Make sure chromosome names in BAM and VCF match (e.g., 'chr1' vs '1')." >&2
    exit 1
fi
CHUNK_SIZE=$((TOTAL_SIZE / NUM_CHUNKS))
echo "Number of threads: ${nprocs}"
echo "Target chunk size: ${CHUNK_SIZE} bp for ${NUM_CHUNKS} total chunks."

# Step 4: Iterate over each chromosome to generate the list of regions.
REGIONS=""
while read -r chrom chrom_len; do
    current_pos=1
    while [[ "$current_pos" -lt "$chrom_len" ]]; do
        end_pos=$((current_pos + CHUNK_SIZE - 1))
        # Ensure the chunk does not extend beyond the end of the chromosome.
        if [[ "$end_pos" -gt "$chrom_len" ]]; then
            end_pos=$chrom_len
        fi
        # Add the region (format chr:start-end) to our list of tasks.
        REGIONS="${REGIONS} ${chrom}:${current_pos}-${end_pos}"
        current_pos=$((end_pos + 1))
    done
done <<< "$GENOME_GEOMETRY"

echo "Generated $(echo "$REGIONS" | wc -w | xargs) tasks for GNU Parallel."

# ==============================================================================
# --- Phase 2: Parallel Execution (Map) ---
# ==============================================================================

# The worker function for each thread, now adapted for a region.
scatter_and_run() {
    local region="$1"
    local region_safe_name
    region_safe_name=$(echo "$region" | tr ':-' '_')
    
    local temp_output="${TMPDIR}/pileup_${region_safe_name}.csv"
    local temp_nbam="${TMPDIR}/normal_${region_safe_name}.bam"
    local temp_tbam="${TMPDIR}/tumor_${region_safe_name}.bam"
    local temp_vcf="${TMPDIR}/snps_${region_safe_name}.vcf.gz"

    # Extract data for THIS SPECIFIC REGION
    "${SAMTOOLS_EXE}" view -b "${normal_bam}" "${region}" -o "${temp_nbam}"
    "${SAMTOOLS_EXE}" view -b "${tumor_bam}" "${region}" -o "${temp_tbam}"
    "${BCFTOOLS_EXE}" view --regions "${region}" "${snp_vcf}" -Oz -o "${temp_vcf}"

    # Index the temporary files
    "${SAMTOOLS_EXE}" index "${temp_nbam}"
    "${SAMTOOLS_EXE}" index "${temp_tbam}"
    "${BCFTOOLS_EXE}" index "${temp_vcf}"

    # Run snp-pileup
    "${SNP_PILEUP_EXE}" --pseudo-snps="${pseudo_snps}" -q "${mapq}" -Q "${baseq}" ${count_orphans} "${temp_vcf}" "${temp_output}" "${temp_nbam}" "${temp_tbam}"
}

export -f scatter_and_run
export normal_bam tumor_bam snp_vcf TMPDIR SNP_PILEUP_EXE SAMTOOLS_EXE BCFTOOLS_EXE pseudo_snps mapq baseq count_orphans

echo "Starting parallel processing with ${nprocs} threads..."
parallel --jobs "${nprocs}" scatter_and_run ::: ${REGIONS}
echo "Parallel processing finished."

# ==============================================================================
# --- Phase 3: Gathering Results (Reduce) ---
# ==============================================================================

echo "Concatenating, sorting, and compressing results..."
FIRST_FILE=$(find "$TMPDIR" -name '*.csv' -print -quit)
if [ -z "$FIRST_FILE" ]; then
    echo "Error: No pileup files were generated." >&2
    exit 1
fi

# A robust pipeline to gather the results
{
    # Print the header from the first file.
    head -n 1 "$FIRST_FILE";
    # Concatenate the content of all other files, skipping their headers.
    tail -q -n +2 "${TMPDIR}"/*.csv;
} | sort -k1,1V -k2,2n | gzip > "$output_pileup"

echo "Script finished successfully. Final output is in ${output_pileup}"
