#!/usr/bin/python
# script to modelize signatures values
# version 1 - 11-12-2012
# Usage signature_sampling.py <bowtie input> <output> <bowtie index>

import sys, subprocess, random
from collections import defaultdict # required for some SmRNAwindow attributes (readDic)
from numpy import mean, std # required for some SmRNAwindow methods
from smRtools import *

def read_collect(line):
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

def report_sig():
  for item in objDic:
    item_reads = objDic[item].readcount(23,29)
    read_density = item_reads / float(objDic[item].size)
    zscore_sig = objDic[item].z_signature(23,29,23,29, range(1,26) )[10]
    prob_sig = objDic[item].hannon_signature(23,29,23,29, range(1,26) )[10] * 100
    print >> OUT, "%s\t%s\t%s\t%s\t%s\t%s" % (n+1, item, item_reads, read_density, zscore_sig, prob_sig)


fasta_dic = get_fasta (sys.argv[3])
F = open (sys.argv[1], "r") # F is the bowtie output taken as input
Fullfile = F.read().splitlines() # a list of lines in RAM
F.close()

OUT = open (sys.argv[2], "w")
print >> OUT, "%s\t%s\t%s\t%s\t%s\t%s" % ("Total_Reads", "item", "item_reads", "read_density", "zscore_sig", "prob_sig")

for round in range(200):
  random.shuffle(Fullfile) # to mix well !
  objDic = {}
  for n, line in enumerate(Fullfile):
    read_collect(line) # feeds objDic
    if (n+1)%1500 == 0:
      report_sig()

OUT.close()
