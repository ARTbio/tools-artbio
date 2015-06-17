#!/usr/bin/python
# script for computing small RNA phasing from a bowtie output
# version 1 4-1-2013
# Usage small_RNA_phasing.py <bowtie input> <minsize query> <maxsize query> <min scope> <max scope> <output> <bowtie index>

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

minquery = int(sys.argv[2])
maxquery = int(sys.argv[3])
minscope = int(sys.argv[4])
maxscope = int(sys.argv[5]) + 1
scope = range (minscope, maxscope)
sizerange = range(minquery, maxquery+1)
general_phasing_table = dict ([(i,0) for i in scope])
######
OUT = open (sys.argv[6], "w")
for x in (objDic):
  local_phasing_table = objDic[x].phasing( sizerange, scope )
  try:
    for phase_offset in local_phasing_table.keys():
      general_phasing_table[phase_offset] = general_phasing_table.get(phase_offset, 0) + local_phasing_table[phase_offset]
  except:
    pass

print >> OUT, "Phase_offset\tprobability"

for phase_offset in sorted(general_phasing_table):
  print >> OUT, "%i\t%f" % (phase_offset, general_phasing_table[phase_offset])

OUT.close()


