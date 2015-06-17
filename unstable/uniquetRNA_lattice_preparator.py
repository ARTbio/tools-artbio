#!/usr/bin/python
# python parser module for lattice preparation from  bowtie 23/6/2012
# specifically for tRNA "canonical" (Anahi's project)
# version 1
# rewritting of the module based on the smRtool.py class
# Usage uniquetRNA_lattice_preparator.py <bowtie_out> <output1 file> <norm_factor>

import sys, subprocess
from collections import defaultdict
from smRtools import *
from numpy import mean, median, std

class LatticeRNA (SmRNAwindow):
  '''overloading of the smRNAwindow class for objects with only forward reads (typically mRNA matched by reads)'''

  def readmap (self):
    readmap = {}
    for offset in self.readDict.keys():
      readmap[offset] = len(self.readDict[offset])
    return readmap

  def normalizedreadmap (self):
    MaxOffset = self.size
    readmap = {}
    thevalues=[]
    for offset in self.readDict.keys():
      thevalues.append(len(self.readDict[offset]))
    try: MaxValue = max(thevalues)
    except: MaxValue = 0
    for offset in self.readDict:
      readmap[offset/float(MaxOffset)] = len(self.readDict[offset])/float(MaxValue)
    return readmap    

  def meansizeatoffset (self, estimator_function, offset):
    return estimator_function(self.readDict[offset])


  def meansizemap (self, estimator_function):
    meansizedic = {}
    for offset in self.readDict.keys():
      meansizedic[offset] = estimator_function(self.readDict[offset])
    return meansizedic

  def medianesizemap (self):
    medianesizedic = {}
    for offset in self.readDict.keys():
      medianesizedic[offset] = median(self.readDict[offset])
    return medianesizedic

  def density (self):
    '''method to output the read coverage by position in the mir'''
    map = [0 for i in range (len(self.sequence))]
    for offset, size in self.dicmap:
      for i in range (offset, offset+size):
        map[i] += self.dicmap[(offset,size)]
    return map

  def normalized_density (self):
    map = self.density ()
    maximum = float (max (map) ) or 1
    length = float (len (map) ) or 1
    Total_NoR = self.mircount()
    output = ["mir\tcoordinate\tdensity\tNoR"]
    for i, D in enumerate (map):
      output.append("%s\t%s\t%s\t%s" % (self.name, (i+1)/length, D/maximum, Total_NoR))
    return "\n".join(output)

Cannonical_tRNA = LatticeRNA("Cannonical tRNA")
   
F = open (sys.argv[1], "r")
for line in F:
  fields = line.split()
  name = "Cannonical tRNA"
  offset= int(fields[2])
  sequence= fields[3]
  Cannonical_tRNA.addread("+", offset, len(sequence))
F.close()

norm_factor = sys.argv[3]
norm_factor = float(norm_factor)
F = open (sys.argv[2], "w")
print >> F, "gene\toffset\tcount\tmedianesize\ttotal_count"
for offset in sorted(Cannonical_tRNA.readDict):
  print >> F, "%s\t%s\t%s\t%s\t%s" % (Cannonical_tRNA.gene, offset, len(Cannonical_tRNA.readDict[offset])*norm_factor, int(Cannonical_tRNA.meansizeatoffset(median, offset)), Cannonical_tRNA.forwardreadcount()*norm_factor )
F.close()

