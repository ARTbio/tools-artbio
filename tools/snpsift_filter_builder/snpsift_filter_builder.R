# Create Filtering Input for variant selection with SnpSift tool based on sample names and sample phenotypes
# Exemple
# nolint start
# "( AF<0.01 | AF =='.' ) & isVariant( GEN[0] ) & isVariant( GEN[2] ) & countRef()>=1"
# nolint end
# Command line
# "Rscript snpsift_filter_builder --description_table_file desc_test.txt --header_desc T --sample_col 1 --phenotype_col 2 --vcf test.vcf --normal normal --patient patient --af 0.01 --count_ref 2"
# Looking for variants in patient samples
# with allele frequency smaller than 0.01
# and a minimal of 1 reference allele in the genotypes

options(show.error.messages = F,
        error = function() {
            cat(geterrmessage(), file = stderr());q("no", 1, F)
            }
        )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
warnings()

library(optparse)

# Arguments
option_list <- list(
  make_option(
    "--description_table_file",
      default = NA,
      type = "character",
      help = "Table describing samples and phenotypes"
    ),
  make_option(
    "--header_desc",
    default = T,
      type = "logical",
      help = "True/False Header present in description table or not"
    ),
  make_option(
    "--sample_col",
    default = 1,
    type = "integer",
    help = "Column number of sample names"
    ),
  make_option(
      "--phenotype_col",
      default = 2,
      type = "integer",
      help = "Column number of sample phenotypes"
  ),
  make_option(
    "--vcf",
    default = NA,
    type = "character",
    help = "VCF file name containing data to filter with sample names"
  ),
  make_option(
    "--normal",
    default = "normal",
    type = "character",
    help = "Phenotype value of normal samples"
  ),
  make_option(
    "--patient",
    default = "patient",
    type = "character",
    help = "Phenotype value of not normal samples"
  ),
  make_option(
    "--af",
    default = NA,
    type = "numeric",
    help = "Maximum AF allele frequency value"
  ),
  make_option(
    "--count_ref",
    default = -1,
    type = "integer",
    help = "Minimum number of reference alleles"
  ),
  make_option(
    c("-o", "--out"),
    default = "filter_line.txt",
    type = "character",
    help = "Output name [default : '%default' ]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

################Manage input data ####################


description_table_file <- opt$description_table_file
header_desc <- opt$header_desc
sample_col <- opt$sample_col
pheno_col <- opt$phenotype_col
vcf_file <- opt$vcf
norm_value <- opt$normal
patient_value <- opt$patient
af <- opt$af
count_ref <- opt$count_ref
output_file <- opt$out


print(c(description_table_file, header_desc, sample_col, pheno_col, vcf_file,
        norm_value, patient_value, af, count_ref, output_file))

# Read Sample description table and get sample names and phenotypes
desc_table <- read.table(file = description_table_file,
                         sep = "\t",
                         header = header_desc)

data_table <- data.frame(sampleID = desc_table[, sample_col],
                         phenotype = desc_table[, pheno_col])

# Check presence of normal and not normal samples in table description
if (length(grep(norm_value, data_table$phenotype)) < 1) warning(
    paste0("There are no normal samples in description table with phenotype ", normal_value))
if (length(grep(patient_value, data_table$phenotype)) < 1) warning(
    paste0("There are no patient samples in description table with phenotype ", patient_value))


# Get VCF #CHROM line containing the names and ordering of samples in VCF
chrom_line <- system(paste0("grep -m 1 ^#CHROM ", vcf_file), intern = T)
samples <- unlist(strsplit(gsub(".*FORMAT\t", "", chrom_line), "\t"))
samples_data <- data.frame(sampleID = samples, position = seq(0, length(samples) - 1, by = 1))

# Verifying that all samples in VCF are in description table
 if (length(intersect(samples_data$sampleID, data_table$sampleID)) != length(samples)) {
    stop(paste(paste0(setdiff(samples_data$sampleID, data_table$sampleID), collapse = ","),
               " sample(s) not present in description table"))
}

# Merging data_table and position in VCF
data_table <- merge(samples_data, data_table)

# Get filter line
filter_line <- NULL

# Add AF
if (! is.na(af)) filter_line <- paste0("( AF<", af, " | AF == '.' )")

# add variants
if (length(grep(patient_value, data_table$phenotype)) > 0) {
    for (i in grep(patient_value, data_table$phenotype)) {
        filter_line <- ifelse(is.null(filter_line),
                              paste0(" isVariant( GEN[", data_table$position[i], "] ) "),
                              paste0(filter_line, " & isVariant( GEN[", data_table$position[i], "] ) "))
  }
}

# add countRef
if (count_ref != -1) {
    filter_line <- ifelse(is.null(filter_line),
                          paste0("countRef()>= ", count_ref),
                          paste0(filter_line, " & countRef()>= ", count_ref))
    } else {
    if (length(grep(norm_value, data_table$phenotype)) > 1)
        filter_line <- ifelse(is.null(filter_line),
                         paste0("countRef()>= ", length(grep(norm_value, data_table$phenotype)) - 1),
                         paste0(filter_line, " & countRef()>= ", length(grep(norm_value, data_table$phenotype)) - 1))
}

#Return filter line
print(filter_line)
write(filter_line, file = output_file)
