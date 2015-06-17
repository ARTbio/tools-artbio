#!/usr/bin/python
# script for outputing pairABLE reads from a bowtie output
# version 1 version 3-12-2012
# Usage pairable.py <bowtie input> <minsize query> <maxsize query> <minsize target> <maxsize target> <output> <bowtie index>

import sys, subprocess
from collections import defaultdict
from smRtools import get_fasta, antipara, RNAtranslate, SmRNAwindow

fasta_dic = get_fasta (sys.argv[7])
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

OUT = open (sys.argv[6], "w")


for x in objDic:
  sequence_list= objDic[x].newpairable_bowtie( 10, sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
  if sequence_list:
    for line in sequence_list:
      print >> OUT, line

OUT.close()
