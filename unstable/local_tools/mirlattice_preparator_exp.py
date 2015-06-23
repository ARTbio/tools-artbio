#!/usr/bin/python
# python parser module for lattice preparation from  bowtie 23/6/2012
# version 2
# rewritting of the module based on the smRtool.py class
# Usage mirlattice_preparator.py <bowtie_out> <index> <output1 file>

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
    try:
      MaxOffset = max(self.readDict)
    except:
      return self.readmap()
    readmap = {}
    thevalues=[]
    for offset in self.readDict.keys():
      thevalues.append(len(self.readDict[offset]))
    MaxValue = max(thevalues)
    for offset in self.readDict:
      readmap[offset/float(MaxOffset)] = len(self.readDict[offset])/float(MaxValue)
    return readmap    

  def meansizeatoffset (self, estimator_function, offset):
    return estimator_function(self.readDict[offset])

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

ItemDic = get_fasta (sys.argv[2])
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

F = open (sys.argv[3], "w")
print >> F, "gene\tnormOffset\tsize"
for item in sorted(ObjectDic):
  for offset in sorted(ObjectDic[item].readDict):
    normoffset = offset / float(ObjectDic[item].size)
    for size in ObjectDic[item].readDict[offset]:
      print >> F, "%s\t%s\t%s" % (item, normoffset, size)
F.close()

