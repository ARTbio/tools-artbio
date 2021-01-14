#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DNAcopy")
#library(DNAcopy)
#install.packages("saasCNV", dependencies = TRUE)
#library(saasCNV)

######################
##### Modified function vcf2txt because MQ is not available in our VCFs and mandatory in "standard analysis" to filter reads 
### (MQ.cutoff>= 30 by default)
### Think about another way to filter variants???

vcf2txt2<-function (vcf.file, normal.col = 10, tumor.col = 11,MQ.cutoff=30) 
{
  inputFile <- file(vcf.file, "r")
  while (length(lines <- readLines(inputFile, n = 1, warn = FALSE)) > 
         0) {
    if (length(grep("\\#CHROM", lines)) == 1) 
      break
  }
  
  
  header <- sub("\\#", "", lines)
  header <- strsplit(header, "\t")[[1]]
  header[c(normal.col, tumor.col)] <- c("Normal", "Tumor")
  vcf <- read.table(inputFile, header = FALSE, sep = "\t", 
                    as.is = TRUE, comment.char = "")
  names(vcf) <- header
  vcf <- vcf[, c(1:9, normal.col, tumor.col)]
  close(inputFile)
  vcf$CHROM <- ifelse(grepl("^chr", vcf$CHROM), vcf$CHROM, paste0("chr", vcf$CHROM))
  chrs <- paste0("chr", c(1:22, "X", "Y"))
  idx <- with(vcf, CHROM %in% chrs & FILTER == "PASS")
  vcf <- vcf[idx, ]
  idx <- grep(",", vcf$ALT)
  if (length(idx) >= 1) 
    vcf <- vcf[-idx, ]
  normal.RD <- strsplit(vcf$Normal, ":")
  normal.RD <- do.call(rbind, normal.RD)
  vcf$Normal.GT <- normal.RD[, 1]
  ad.col=4
  normal.AD <- do.call(rbind, strsplit(normal.RD[, ad.col], ","))
  vcf$Normal.REF.DP <- as.integer(normal.AD[, 1])
  vcf$Normal.ALT.DP <- as.integer(normal.AD[, 2])
  tumor.RD <- strsplit(vcf$Tumor, ":")
  tumor.RD <- do.call(rbind, tumor.RD)
  vcf$Tumor.GT <- tumor.RD[, 1]
  tumor.AD <- do.call(rbind, strsplit(tumor.RD[, ad.col], ","))
  vcf$Tumor.REF.DP <- as.integer(tumor.AD[, 1])
  vcf$Tumor.ALT.DP <- as.integer(tumor.AD[, 2])
  vcf <- vcf[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                 "INFO", "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", 
                 "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]
  idx <- which(vcf$Normal.GT == "./." | vcf$Tumor.GT == "./.")
  if (length(idx) >= 1) 
    vcf <- vcf[-idx, ]
  info <- strsplit(vcf$INFO, ";")
  extract.MQ <- function(x) {
    idx.MQ <- grep("^MQ\\=", x)
    if (length(idx.MQ) == 1) 
      return(x[idx.MQ])
    else return(NA)
  }
  info1 <- sapply(info, FUN = extract.MQ)
  
  if(!is.na(sum(as.numeric(sub("^MQ\\=", "", info1))))){
    vcf$MQ <- as.numeric(sub("^MQ\\=", "", info1))
    
    vcf <- vcf[vcf$MQ >= MQ.cutoff, ]
    vcf <- vcf[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                   "MQ", "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", 
                   "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]
  }
  a <- vcf[, c("CHROM", "POS")]
  idx <- which(duplicated(a))
  if (length(idx) >= 1) 
    vcf <- vcf[-idx, ]
  return(vcf)
}



############
### Fonction qui rÃ©alise l'analyse CNV pour 1 Ã©chantillon
### input nÃ©cessaire
## hist_id_vcf : numÃ©ro du vcf dans l'historique
## hist_id_refGene : numÃ©ro du refGene.txt avec header dans l'historique
## normal.col/tumor.col : numÃ©ro de colonne qui contient les infos de l'Ã©chantillon normal/tumor dans le vcf
## output :
## 4 figures au format png et 1 tableau au format tabular qui contient les CNV associÃ©s Ã  leur pvalue 

sampleCNV <- function( hist_id_vcf = hist_id_vcf, normal.col=normal.col, tumor.col=tumor.col, sample.name=sample.name , hist_id_refGene = hist_id_refGene ){
  
  vcf.file<-gx_get(hist_id_vcf)
  
  # Convert vcf to table and filter "PASS" variants
  vcf_table <- vcf2txt2(vcf.file=vcf.file, normal.col=normal.col, tumor.col=tumor.col)
  
  
  # Transform read depth information into log2ratio and log2mBAF
  seq.data <- cnv.data(vcf=vcf_table, min.chr.probe=100, verbose=TRUE)
  
  # Perform joint segmentation on log2ratio and log2mBAF dimensions. The function joint.segmentation
  # outputs the starting and ending points of each CNV segment as well as some summary statistics
  seq.segs <- joint.segmentation(data=seq.data, min.snps=10, global.pval.cutoff=1e-4, max.chpts=30, verbose=TRUE)
  
  
  # merge adjacent segments, for which the median values in either or both dimensions are not substantially different
  
  #seq.segs.merge <- merging.segments(data=seq.data, segs.stat=seq.segs, use.null.data=TRUE,N=1000, maxL=20000, merge.pvalue.cutoff=0.05, verbose=TRUE)
  seq.segs.merge <- merging.segments(data=seq.data, segs.stat=seq.segs, use.null.data=TRUE, merge.pvalue.cutoff=0.05, verbose=TRUE)
  
  # # Plot result for one chromosome
  # chrom=7
  # png(paste0("diag_",sample.name,"_chr",chrom,".png"),width=1400,height=700 )
  # diagnosis.seg.plot.chr(data=seq.data, segs=seq.segs,sample.id=sample.name,chr=7, cex=0.3)
  # dev.off()
  # gx_put(filename = paste0("diag_",sample.name,"_chr",chrom,".png"),file_type = "png")
  
  lapply(c(1:22,"X","Y"),function(x){
    chrom=x
    png(paste0("diag_",sample.name,"_merged_chr",chrom,".png"),width=1400,height=700 )
    diagnosis.seg.plot.chr(data=seq.data, segs=seq.segs.merge,sample.id=paste(sample.name,"merged segments"),chr=chrom, cex=0.3)
    dev.off()
    gx_put(filename = paste0("diag_",sample.name,"_merged_chr",chrom,".png"),file_type = "png")
  })
  
  ### CNV Analysis on genome
  
  seq.cnv <- cnv.call(data=seq.data, sample.id=sample.name, segs.stat=seq.segs.merge, pvalue.cutoff=0.05)
  
  gene.anno.file=gx_get(hist_id_refGene)
  gene.anno <- read.delim(file=gene.anno.file, as.is=TRUE, comment.char="")
  seq.cnv.anno <- reannotate.CNV.res(res=seq.cnv, gene=gene.anno, only.CNV=TRUE)
  write.table(seq.cnv.anno, file = paste0("CNV_",sample.name,".txt"), sep="\t",quote=F,row.names=F)
  gx_put(filename = paste0("CNV_",sample.name,".txt"),file_type = "tabular")
  ## Plot results on genome
  
  png(paste0("genomeCNV_",sample.name,".png"),width=1400,height=700 )
  genome.wide.plot(data=seq.data, segs=seq.cnv, sample.id=sample.name, chrs=sub("^chr","",unique(seq.cnv$chr)),cex=0.3)
  dev.off()
  gx_put(filename = paste0("genomeCNV_",sample.name,".png"),file_type = "png")
  
  png(paste0("diagCNV_",sample.name,".png"),width=700,height=700 )
  diagnosis.cluster.plot(segs=seq.cnv, chrs=sub("^chr","",unique(seq.cnv$chr)), min.snps=10, max.cex=3, ref.num.probe=1000)
  dev.off()
  gx_put(filename = paste0("diagCNV_",sample.name,".png"),file_type = "png")
  
  #gx_save()
  
}



####### Lancement de l'analyse avec une boucle for (je sais!!)
## lignes de commandes pour lancer une analyse sur 1 collection de vcf
## input :

## hist_id_tablePos : numÃ©ro du tableau dans l'historique
## 		tableau de 2 colonnes dans l'historique galaxy
##			colonne 1: nom de l'Ã©chantillon
##			colonne 2: position de l'Ã©chantillon dans la collection

## hist_id_vcf_collection : numÃ©ro de la collection de vcf dans l'historique
## hist_id_refGene : numÃ©ro du refGene.txt avec header dans l'historique
## normal.col/tumor.col : numÃ©ro de colonne qui contient les infos de l'Ã©chantillon normal/tumor dans le vcf
## 

hist_id_tablePos <- 38
hist_id_refGene <- 79
hist_id_vcf_collection <- 31
normal.col = 9+1
tumor.col = 9+2




### Code Ã  lancer dans Rstudio Ã  partir de l'historique : https://usegalaxy.sorbonne-universite.fr/u/fanconi/h/cnv-analysis 
tabSamplesIDs <- read.delim(gx_get(hist_id_tablePos))

nbSamples <- nrow(tabSamplesIDs)

# for(i in 1:nbSamples){
#   print(paste0("CNV analysis for ", tabSamplesIDs$sample[i]))
#   sampleCNV(hist_id_vcf = hist_id_vcf_collection-nbSamples+tabSamplesIDs$collection_position[i]-1, 
#             normal.col=normal.col, tumor.col=tumor.col, 
#             sample.name=tabSamplesIDs$sampleID[i], 
#             hist_id_refGene = hist_id_refGene )
#   
# }
# 
lapply(3:nbSamples, function(i){
  print(paste0("CNV analysis for ", tabSamplesIDs$sample[i]))
  sampleCNV(hist_id_vcf = hist_id_vcf_collection-nbSamples+tabSamplesIDs$collection_position[i]-1,
            normal.col = normal.col, tumor.col = tumor.col,
            sample.name = tabSamplesIDs$sampleID[i],
            hist_id_refGene = hist_id_refGene )
})

