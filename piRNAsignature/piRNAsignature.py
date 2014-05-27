#!/usr/bin/python
# script for computing overlap signatures from a bowtie output
# Christophe Antoniewski <drosofff@gmail.com>
# Usage piRNAsignature.py <1:input> <2:format of input> <3:minsize query> <4:maxsize query> <5:minsize target> <6:maxsize target> <7:minscope> <8:maxscope> <9:output> <10:bowtie index> <11:procedure option>

import sys, subprocess
from smRtools import *
from collections import defaultdict # test whether it is required

if sys.argv[11] == "--extract_index":
  if sys.argv[2] == "tabular":
    Genome = HandleSmRNAwindows (sys.argv[1],"tabular",sys.argv[10],"bowtieIndex")
  elif sys.argv[2] == "sam":
    Genome = HandleSmRNAwindows (sys.argv[1],"sam",sys.argv[10],"bowtieIndex")
  else:
    Genome = HandleSmRNAwindows (sys.argv[1],"bam",sys.argv[10],"bowtieIndex")
else:
  if sys.argv[2] == "tabular":
    Genome = HandleSmRNAwindows (sys.argv[1],"tabular",sys.argv[10],"fastaSource") 
  elif sys.argv[2] == "sam":
    Genome = HandleSmRNAwindows (sys.argv[1],"sam",sys.argv[10],"fastaSource")
  else:
    Genome = HandleSmRNAwindows (sys.argv[1],"bam",sys.argv[10],"fastaSource")
# this decisional tree may be simplified if sam and bam inputs are treated the same way by pysam

# replace objDic by Genome.instanceDict or... objDic = Genome.instanceDict
objDic = Genome.instanceDict

minquery = int(sys.argv[3])
maxquery = int(sys.argv[4])
mintarget = int(sys.argv[5])
maxtarget = int(sys.argv[6])
minscope = int(sys.argv[7])
maxscope = int(sys.argv[8]) + 1
general_frequency_table = dict ([(i,0) for i in range(minscope,maxscope)])
general_percent_table = dict ([(i,0) for i in range(minscope,maxscope)])


###### for normalized summing of local_percent_table(s)
readcount_dic = {}
Total_read_in_objDic = 0
for item in objDic:
  readcount_dic[item] = objDic[item].readcount(minquery, maxquery)
  Total_read_in_objDic += readcount_dic[item]
######


OUT = open (sys.argv[9], "w")
for x in (objDic):
  local_frequency_table = objDic[x].signature( minquery, maxquery, mintarget, maxtarget, range(minscope,maxscope) )
  local_percent_table = objDic[x].hannon_signature( minquery, maxquery, mintarget, maxtarget, range(minscope,maxscope) )
  try:
    for overlap in local_frequency_table.keys():
      general_frequency_table[overlap] = general_frequency_table.get(overlap, 0) + local_frequency_table[overlap]
  except:
    pass
  try:
    for overlap in local_percent_table.keys():
      general_percent_table[overlap] = general_percent_table.get(overlap, 0) + (1./Total_read_in_objDic*readcount_dic[x]*local_percent_table[overlap])
  except:
    pass

print >> OUT, "overlap\tnum of pairs\tprobability"

# if minquery==mintarget and maxquery==maxtarget: # the division / 2 decision is now directly taken in the method of smRtools (26/11/2013)
for classe in sorted(general_frequency_table):
  print >> OUT, "%i\t%i\t%f" % (classe, general_frequency_table[classe], general_percent_table[classe])
#else:
#  for classe in sorted(general_frequency_table):
#    print >> OUT, "%i\t%i\t%f" % (classe, general_frequency_table[classe], general_percent_table[classe])

OUT.close()
