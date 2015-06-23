#!/usr/bin/python
# script extract a particular genome region and compute a z-score piRNA signature
# version 1 - 22-06-2012
# Usage bowtie_window_analysis.py <bowtie input> <geneID> <Upstream_coordinate> <Downstream_coordinate> <bowtie index> <output>

import sys, subprocess
from collections import defaultdict # required for some SmRNAwindow attributes (readDic)
from numpy import mean, std # required for some SmRNAwindow methods
from smRtools import *

geneID = sys.argv[2]
Upstream_coordinate = int(sys.argv[3])
Downstream_coordinate = int(sys.argv[4])
fasta_dic = get_fasta (sys.argv[5])
geneSequence = fasta_dic[sys.argv[2]][Upstream_coordinate:Downstream_coordinate]
geneObject= SmRNAwindow(geneID, geneSequence)

F = open (sys.argv[1], "r") # F is the bowtie output taken as input
counter = 0
for line in F:
  fields = line.split()
  if  fields[2] != geneID : continue
  polarity = fields[1]
  coordinate = int(fields[3])
  if (coordinate < Upstream_coordinate or coordinate > Downstream_coordinate) : continue
  size = len(fields[4])
  geneObject.addread (polarity, coordinate, size)

F.close()

OUT = open (sys.argv[6], "w")
pipi_z = geneObject.z_signature(23,28,23,28, range(1,26) )
print >> OUT, "pipi signature"
print >> OUT, pipi_z
print >> OUT, "sisi signature"
print >> OUT, geneObject.z_signature(20,22,20,22, range(1,26) )
print >> OUT, "total read analyzed"
print >> OUT, geneObject.readcount()
print >> OUT, "size distribution of these reads"
print >> OUT, geneObject.readsizes()
OUT.close()
