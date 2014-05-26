#!/usr/bin/python
# test script for class instance extraction
# version 1 3-5-2012 uses SmRNAtools class
# Usage classextract_test.py <bowtie input>  <bowtie index> <output>

import sys, subprocess
from smRtools import get_fasta, antipara, RNAtranslate, SmRNAwindow
from collections import defaultdict

geneDic = get_fasta (sys.argv[2])
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
    objDic[gene] = SmRNAwindow(gene, geneDic[gene])
    objDic[gene].addread (polarity, offset, size)
F.close()


OUT = open (sys.argv[3], "w")
width = 10000
step = 10000
for pointer in range(width/2,objDic['2L'].size-(width/2), step):
  print >>OUT,  "%s\t%s" % ("2L", objDic['2L'].sub_z_signature (pointer, width, 24,28,24,28, range(1,26) ))
OUT.close()


