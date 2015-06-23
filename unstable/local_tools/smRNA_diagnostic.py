#!/usr/bin/python
# script to use all the smRtools methods and functions
# version 1 - 29-05-2012
# Usage smRNA_diagnostic.py <bowtie input> <output> <bowtie index>

import sys, subprocess
from collections import defaultdict # required for some SmRNAwindow attributes (readDic)
from numpy import mean, std # required for some SmRNAwindow methods
from smRtools import *

fasta_dic = get_fasta (sys.argv[3])
objDic = {}
F = open (sys.argv[1], "r") # F is the bowtie output taken as input

for line in F:
  fields = line.split()
  polarity = fields[1]
  gene = fields[2]
  offset = int(fields[3])
  size = len (fields[4])
  try:
    objDic[gene].addread (polarity, offset, size)
  except KeyError:
    objDic[gene] = SmRNAwindow(gene, fasta_dic[gene])
    objDic[gene].addread (polarity, offset, size)
F.close()

OUT = open (sys.argv[2], "w")

print >> OUT, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("Id", "total piRNAs (24-29nt)", "forward piRNA ratio", "piRNA density", "z-score piRNA", "piRNA prob (%)", "Ufreq (forward | reverse)", "total siRNAs (20-21nt)", "forward siRNA ratio", "siRNA density", "z-score siRNA", "siRNA prob (%)", "piRNA/siRNA ratio")
for item in objDic:
  pitotal = objDic[item].readcount(24,29)
  try:
    piratio = objDic[item].forwardreadcount(24,29) / float(pitotal)
  except ZeroDivisionError:
    piratio = "NA"
  pidensity = pitotal / float(objDic[item].size)
  pipi_z = objDic[item].signature(24,29,24,29, range(1,26), zscore="yes" )[10]
  pipi_prob = objDic[item].hannon_signature(24,29,24,29, range(1,26) )[10] * 100
  Ufreq = objDic[item].Ufreq_stranded (range(24,30) ) # be carreful to the range edges...
  sitotal = objDic[item].readcount(20,21)
  try:
    siratio = objDic[item].forwardreadcount(20,21) / float(sitotal)
  except ZeroDivisionError:
    siratio = "NA"
  sidensity = sitotal / float(objDic[item].size)
  sipi_z = objDic[item].signature(23,29,20,21, range(1,26), zscore="yes" )[10]
  pisi_prob = objDic[item].hannon_signature(23,29,20,21, range(1,26) )[10] * 100
  try:
    pisi_ratio = float(pitotal)/(pitotal+sitotal)
  except ZeroDivisionError:
    pisi_ratio = "NA"
  print >> OUT, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (item, pitotal, piratio, pidensity, pipi_z, pipi_prob, Ufreq, sitotal, siratio, sidensity, sipi_z, pisi_prob, pisi_ratio )

OUT.close()
