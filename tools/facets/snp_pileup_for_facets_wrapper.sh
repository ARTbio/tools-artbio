#!/bin/bash

# ==============================================================================
#
# snp_pileup_for_facets_wrapper.sh (v3.3 - Final Robust Parallelism)
#
# Description:
#   A robust wrapper script for running snp-pileup in parallel. It works around
#   snp-pileup's inability to filter by region by creating temporary per-chromosome
#   input files for each parallel job.
#
# Author:
#   drosofff (ARTbio) & Gemini (Google)
#
# ==============================================================================

# --- Help Function ---
usage() {
    cat << EOF
Usage: $0 -n <normal.bam> -t <tumor.bam> -v <snps.vcf.gz> -o <output.csv.gz> [OPTIONS]

This script runs snp-pileup in parallel across all chromosomes.

REQUIRED ARGUMENTS:
  -n    Path to the normal BAM file (must be indexed).
  -t    Path to the tumor BAM file (must be indexed).
  -v    Path to the SNP reference VCF file (must be bgzip compressed and tabix indexed).
  -o    Path for the final output pileup file (will be bgzip compressed).

OPTIONS:
  -N    Number of parallel processes to use (default: 1).
  -q    Minimum mapping quality for reads (default: 15).
  -Q    Minimum base quality for bases (default: 20).
  -A    Include anomalous read pairs (flag, default: not set).
  -h    Display this help message and exit.
EOF
    exit 0
}

# --- Robustness settings ---
set -e -u -o pipefail

# --- Argument parsing ---
# (This part is correct and remains unchanged)
normal_bam=""
tumor_bam=""
snp_vcf=""
output_pileup=""
nprocs=1
mapq=15
baseq=20
count_orphans=""

while getopts ":hn:t:v:o:N:q:Q:A" opt; do
    case ${opt} in
        h ) usage ;;
        n ) normal_bam="$OPTARG" ;;
        t ) tumor_bam="$OPTARG" ;;
        v ) snp_vcf="$OPTARG" ;;
        o ) output_pileup="$OPTARG" ;;
        N ) nprocs="$OPTARG" ;;
        q ) mapq="$OPTARG" ;;
        Q ) baseq="$OPTARG" ;;
        A ) count_orphans="-A" ;;
        \? ) echo "Invalid option: -$OPTARG" >&2; usage ;;
        : ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# --- Main logic ---

# Find full paths to executables
SNP_PILEUP_EXE=$(command -v snp-pileup)
SAMTOOLS_EXE=$(command -v samtools)
BCFTOOLS_EXE=$(command -v bcftools)

for exe in SNP_PILEUP_EXE SAMTOOLS_EXE BCFTOOLS_EXE; do
    if [ -z "${!exe}" ]; then
        echo "Error: '${exe%%_*}' executable not found in PATH." >&2
        exit 1
    fi
done
echo "Found snp-pileup executable at: ${SNP_PILEUP_EXE}"

echo "Starting SNP pileup process with ${nprocs} parallel jobs..."

TMPDIR=$(mktemp -d)
trap "rm -rf '${TMPDIR}'" EXIT
echo "Temporary directory created at: ${TMPDIR}"

CHROMS=$("${SAMTOOLS_EXE}" view -H "${normal_bam}" | grep "^@SQ" | cut -f 2 | sed 's/SN://' | grep -Fwf - <(zcat "${snp_vcf}" | grep -v "^#" | cut -f 1 | sort -u) )
echo "Found the following chromosomes to process in both BAM and VCF:"
echo "${CHROMS}"

# --- CORRECT PARALLEL EXECUTION LOGIC ---

# Define the function that will be run in parallel for each chromosome.
# It will inherit the exported variables from the main script.
scatter_and_run() {
    local chrom=$1
    echo "Scattering data for chromosome ${chrom}..."

    local temp_vcf="${TMPDIR}/${chrom}.vcf.gz"
    local temp_nbam="${TMPDIR}/${chrom}.normal.bam"
    local temp_tbam="${TMPDIR}/${chrom}.tumor.bam"
    local temp_output="${TMPDIR}/${chrom}.csv"

    "${BCFTOOLS_EXE}" view -r "${chrom}" "${snp_vcf}" -Oz -o "${temp_vcf}"
    "${SAMTOOLS_EXE}" view -b "${normal_bam}" "${chrom}" > "${temp_nbam}"
    "${SAMTOOLS_EXE}" view -b "${tumor_bam}" "${chrom}" > "${temp_tbam}"

    "${SAMTOOLS_EXE}" index "${temp_nbam}"
    "${SAMTOOLS_EXE}" index "${temp_tbam}"

    echo "Running snp-pileup on chromosome ${chrom}..."
    "${SNP_PILEUP_EXE}" --pseudo-snps=300 -q "${mapq}" -Q "${baseq}" ${count_orphans} "${temp_vcf}" "${temp_output}" "${temp_nbam}" "${temp_tbam}"
}

# Export all necessary variables AND the function so they are available to the sub-shells created by GNU Parallel.
export -f scatter_and_run
export snp_vcf normal_bam tumor_bam TMPDIR SNP_PILEUP_EXE SAMTOOLS_EXE BCFTOOLS_EXE mapq baseq count_orphans

# Run the function in parallel, feeding it the list of chromosomes.
parallel --jobs "${nprocs}" scatter_and_run ::: ${CHROMS}

echo "Parallel processing finished. Concatenating results..."

# gather job outputs and sort chromosome remains the same
FIRST_FILE=$(ls -1v "${TMPDIR}"/*.csv 2>/dev/null | head -n 1)
if [ -z "${FIRST_FILE}" ]; then
    echo "Error: No pileup files were generated." >&2
    exit 1
fi

# Use command grouping { ...; } to pipe the combined output of head and tail.
# This entire pipeline writes the final, sorted, compressed file.
{
    # 1. Print the header once.
    head -n 1 "${FIRST_FILE}";

    # 2. Concatenate all files (skipping their headers).
    tail -q -n +2 "${TMPDIR}"/*.csv;

} | sort -t, -k1,1V -k2,2n | bgzip > "${output_pileup}"
echo "Concatenation, sorting and compression complete."
echo "Final output is in ${output_pileup}"
echo "Script finished successfully. The temporary directory will be removed by the trap."
