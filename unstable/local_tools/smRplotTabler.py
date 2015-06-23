#!/usr/bin/python
# script for dataframe preparation from a bowtie output, to plot in lattice
# version 1 14-7-2012 use of a 8th argument to either get the fasta from a bowtie index or directly from a fasta of the history
# Usage smRplot.py <bowtie input> <minsize query> <maxsize query> <output> <bowtie index> <option tag>

import sys, subprocess
from smRtools import get_fasta, get_fasta_from_history, SmRNAwindow
from collections import defaultdict

if sys.argv[-1] == "--extract_index":
  geneDic = get_fasta (sys.argv[5])
else:
  geneDic = get_fasta_from_history (sys.argv[5])

objDic = {}
OUT = open (sys.argv[4], "w")
F = open (sys.argv[1], "r") # F is the bowtie output taken as input
for line in F:
  fields = line.split()
  polarity = fields[1]
  gene = fields[2]
  offset = int(fields[3])
  size = len (fields[4])
  if not (int(sys.argv[2]) <= size <= int(sys.argv[3]) ):
    continue
  try:
    objDic[gene].addread (polarity, offset, size)
  except KeyError:
    objDic[gene] = SmRNAwindow(gene, geneDic[gene])
    objDic[gene].addread (polarity, offset, size)
F.close()

print >> OUT, "gene\tcoord\tcount\tpolarity\titemcount"
for x in (objDic):
  plottable = objDic[x].readplot()
  itemcount = objDic[x].readcount()
  for line in plottable:
    print >> OUT, "%s\t%s" % (line, itemcount)
OUT.close()


