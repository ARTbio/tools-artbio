#!/usr/bin/python
# script to use all the smRtools methods and functions
# version 1 - 26-02-2013
# Usage correlation_mapper.py <bowtie input> <output> <bowtie index> <reference_window_file_path>

import sys, subprocess
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

for transposon in objDic:
  theresults=objDic[transposon].correlation_mapper (sys.argv[4], 588)
  if theresults:
    print >> OUT, transposon, "Pearson Correlation Hits:"
    for i in theresults:
      print >> OUT, i
    OUT.flush()
OUT.close()
