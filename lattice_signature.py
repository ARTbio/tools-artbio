#!/usr/bin/python
# script for computing overlap signatures from a bowtie output IN LATTICE (no summing process)
# 26-12-2012 drosofff@gmail.com
# Use of a 8th argument to either get the fasta from a bowtie index or directly from a fasta of the history
# Usage lattice_signature.py <bowtie input> <minsize query> <maxsize query> <minsize target> <maxsize target> <output> <bowtie index>

import sys, subprocess
from smRtools import get_fasta, get_fasta_from_history, antipara, RNAtranslate, SmRNAwindow
from collections import defaultdict

if sys.argv[-1] == "--extract_index":
  geneDic = get_fasta (sys.argv[7])
else:
  geneDic = get_fasta_from_history (sys.argv[7])

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

general_frequency_table = dict ([(i,0) for i in range(1,26)])
general_percent_table = dict ([(i,0) for i in range(1,26)])
minquery = int(sys.argv[2])
maxquery = int(sys.argv[3])
mintarget = int(sys.argv[4])
maxtarget = int(sys.argv[5])

OUT = open (sys.argv[6], "w")
print >> OUT, "overlap\tnum of pairs\tprobability\titem"
for x in (objDic):
  local_frequency_table = objDic[x].signature( minquery, maxquery, mintarget, maxtarget, range(1,26) )
  local_percent_table = objDic[x].hannon_signature( minquery, maxquery, mintarget, maxtarget, range(1,26) )
#  if minquery==mintarget and maxquery==maxtarget: # the division / 2 decision is now taken in the method of smRtools (26/11/2013)
  for classe in range(1,26):
    print >> OUT, "%i\t%i\t%f\t%s" % (classe, local_frequency_table[classe], local_percent_table[classe], x)
#  else:
#    for classe in range(1,26):
#      print >> OUT, "%i\t%i\t%f\t%s" % (classe, local_frequency_table[classe], local_percent_table[classe], x)
OUT.close()


