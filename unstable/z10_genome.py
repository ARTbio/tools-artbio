#!/usr/bin/python
# script for computing z10 signature from a bowtie genomic output
# version 1 8-5-2012 using SmRNAtools class imported from smRtools
# Usage z10_genome.py <bowtie input> <windowsize> <outputpipi> <outputpisi> <outputsisi> <bowtie index>

import sys, subprocess
from smRtools import get_fasta, antipara, RNAtranslate, SmRNAwindow
from collections import defaultdict
from numpy import mean, std

def split_len(seq, length):
  return [seq[i:i+length] for i in range(0, len(seq), length)]


geneDic = get_fasta (sys.argv[6])
objDic = defaultdict(dict)
windowsize = int(sys.argv[2])

for chrom in geneDic:
  for windowindex, sequencewindow in enumerate (split_len(geneDic[chrom], windowsize)):
    objDic[chrom][windowindex * windowsize] = SmRNAwindow(chrom, sequencewindow, windowindex * windowsize)

F = open (sys.argv[1], "r") # F is the bowtie output taken as input
for line in F:
  fields = line.split()
  polarity = fields[1]
  chrom = fields[2]
  offset = int(fields[3])
  size = len (fields[4])
  objDic[chrom][offset/windowsize*windowsize].addread (polarity, offset, size)
F.close()

PIPI = open (sys.argv[3], "w")
PISI = open (sys.argv[4], "w")
SISI = open (sys.argv[5], "w")

print >> PIPI, 'track type=wiggle_0 name="Z-score 10nt overlap" description="pipi z10_signature" visibility=full color=0,0,255 yLineMark=1 yLineOnOff=on  priority=4'
print >> PISI, 'track type=wiggle_0 name="Z-score 10nt overlap" description="pisi z10_signature" visibility=full color=0,0,255 yLineMark=1 yLineOnOff=on  priority=4'
print >> SISI, 'track type=wiggle_0 name="Z-score 10nt overlap" description="sisi z10_signature" visibility=full color=0,0,255 yLineMark=1 yLineOnOff=on  priority=4'

for chrom in objDic:
  print >>PIPI, "variableStep chrom=chr%s span=%s" % (chrom, windowsize)
  print >>PISI, "variableStep chrom=chr%s span=%s" % (chrom, windowsize)
  print >>SISI, "variableStep chrom=chr%s span=%s" % (chrom, windowsize)
  for windowoffset in sorted(objDic[chrom]):
    value = objDic[chrom][windowoffset].z_signature(24,28,24,28, range(1,26) )[10]
    if value > 2:
      print >>PIPI, "%s\t%.1f" % (windowoffset+1, value)
    value = objDic[chrom][windowoffset].z_signature(20,22,24,28, range(1,26) )[10]
    if value > 2:
      print >>PISI, "%s\t%.1f" % (windowoffset+1, value)
    value = objDic[chrom][windowoffset].z_signature(20,22,20,22, range(1,26))[10]
    if value > 2:
      print >>SISI, "%s\t%.1f" % (windowoffset+1, value)

PIPI.close()
PISI.close()
SISI.close()


