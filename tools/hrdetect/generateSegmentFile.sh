######   Analyse HRDetect - préparation des fichiers input
### Basée sur https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html

## En 3 étapes
# Step1 unix command line : preprocessing des data avec sequenza-utils
# Step2 R code : préparation du fichier de segmentation
# step3 R script : finalisation du fichier segments.tsv avec annotation des lohtype  

## Step1 preprocessing des data 
# Installation de sequenza-utils

pip install sequenza-utils


# Processing du fasta de référence (à faire 1 fois pour hg19 et pour hg38)

sequenza−utils gc_wiggle  --fasta hg19.fa -o hg19.gc50Base.wig.gz −w 50

# Production du segmentation file seqz à partir d'un VCF 
# Ici on prend un output de varscan contenant seulement les SNPs (pas indel) : 
# https://usegalaxy.sorbonne-universite.fr/u/fanconi/h/vcfs-for-hrdetect-input-files-preparation

sequenza-utils snp2seqz -v sample.vcf  -gc hg19.gc50Base.wig.gz -o sample.seqz.gz 

# Post-process by binning the original seqz file:
 
sequenza-utils seqz_binning --seqz sample.seqz.gz -w 50 -o sample.binned.seqz.gz
  
 
 #### Step2 préparation du fichier de segmentation
 
 Rscript segmentation_sequenza.R -i sample.binned.seqz.gz -s sample -d test
 
 ## Step3 
Rscript sequenza_to_hrdtools_input.R -i test/sample_segments.txt -s test/sample_alternative_solutions.txt -o hdrdetect_segment.tsv

 

