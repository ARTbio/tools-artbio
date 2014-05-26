#!/usr/bin/python
# python parser module for lattice preparation from  bowtie 23/6/2012
# version 3 16-4-2014
# Usage mirlattice_preparator.py <bowtie_out> <output file> <norm_factor> <bowtie index> <option tag>

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

if sys.argv[-1] == "--extract_index":
  ItemDic = get_fasta (sys.argv[-2])
else:
  ItemDic = get_fasta_from_history (sys.argv[-2])

ObjectDic = {}
for item in ItemDic:
  ObjectDic[item] = LatticeRNA(item, ItemDic[item])
   
F = open (sys.argv[1], "r")
for line in F:
  fields = line.split()
  name = fields[1]
  offset= int(fields[2])
  sequence= fields[3]
  ObjectDic[name].addread("+", offset, len(sequence))
F.close()
norm_factor = sys.argv[3]
norm_factor = float(norm_factor)
F = open (sys.argv[2], "w")
print >> F, "gene\toffset\tcount\tnormOffset\tnormCount\tmedianesize\ttotal_count"
for item in sorted(ObjectDic):
  for offset, normoffset in zip (sorted(ObjectDic[item].readDict), sorted(ObjectDic[item].normalizedreadmap()) ):
    print >> F, "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (item, offset, len(ObjectDic[item].readDict[offset])*norm_factor, normoffset, ObjectDic[item].normalizedreadmap()[normoffset], int(ObjectDic[item].meansizeatoffset(median, offset)), ObjectDic[item].forwardreadcount()*norm_factor )
F.close()

