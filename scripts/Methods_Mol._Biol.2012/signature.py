#!/usr/bin/env python
# computes overlap signatures from a bowtie or SAM input
# version Met Mol Biol 06-2-2012
# Usage signature.py <bowtie or sam input> <minsize> <maxsize> <minscope> <maxscope> <output>

import sys
from collections import defaultdict
from numpy import mean, std

class SmRNAwindow:
  def __init__(self, gene):
    self.readDict = defaultdict(list) # "dictionary of lists" structure {+/-offset:[size1, size2, ...], ...}
  def addread (self, polarity, offset, size):
    if polarity == "+":
      self.readDict[offset].append(size)
    else:
      self.readDict[-(offset + size -1)].append(size)
    return
  def readcount (self, lim_inf=0, lim_sup=1000):
    n=0
    for offset in self.readDict:
      licenced = [i for i in self.readDict[offset] if (i>=lim_inf and i<= lim_sup)]
      n+=len(licenced)
    return n
  def count_pairs (self, minsize, maxsize, scope):
    size_range = range (minsize, maxsize+1)
    Query_table = {}
    Target_table = {}
    frequency_table = dict ([(i, 0) for i in scope])
    for offset in self.readDict:
      for size in self.readDict[offset]:
        if size in size_range:
          Query_table[offset] = Query_table.get(offset, 0) + 1
          Target_table[offset] = Target_table.get(offset, 0) + 1
    for offset in Query_table:
      for i in scope:
        frequency_table[i] += min(Query_table[offset], Target_table.get(-offset -i +1, 0))
    for i in frequency_table:
      frequency_table[i] = frequency_table[i] / 2 # the number of smRNA PAIRS, not of the paired smRNA.
    return frequency_table
  def overlap_probability (self, minsize, maxsize, scope):
    size_range = range (minsize, maxsize+1)
    Query_table = {}
    Target_table = {}
    Total_Query_Numb = 0
    general_frequency_table = dict ([(i,0) for i in scope])
    for offset in self.readDict:
      for size in self.readDict[offset]:
        if size in size_range:
          Query_table[offset] = Query_table.get(offset, 0) + 1
          Target_table[offset] = Target_table.get(offset, 0) + 1
          Total_Query_Numb += 1
    for offset in Query_table:
      frequency_table = dict ([(i,0) for i in scope])
      number_of_targets = 0
      for i in scope:
        frequency_table[i] += Query_table[offset] *  Target_table.get(-offset -i +1, 0)
        number_of_targets += Target_table.get(-offset -i +1, 0)
      for i in scope:
        try:
          general_frequency_table[i] += (1. / number_of_targets / Total_Query_Numb) * frequency_table[i]
        except ZeroDivisionError :
          continue
    return general_frequency_table      

def load_input (input_file):
  F=open(input_file)
  sampleline = F.readline()
  F.close()
  samplefields=sampleline.split()
  if len(samplefields) < 2 : #alignment format are tabulated with at least two fields. 
    print "error: invalid input format"
    sys.exit()
  if samplefields[1] in ["+", "-"]: # standard tabular bowtie format detected
    F = open (input_file, "r")
    for line in F:
      fields = line.split()
      polarity = fields[1]
      gene = fields[2]
      offset = int(fields[3]) + 1 # to shift on 1-based coordinates
      size = len (fields[4])
      try:
        objDic[gene].addread (polarity, offset, size)
      except KeyError:
        objDic[gene] = SmRNAwindow(gene)
        objDic[gene].addread (polarity, offset, size)
    F.close()
  elif samplefields[0][0] == "@": # SAM format detected
    F = open (input_file, "r")
    for line in F:
      if line[0] == "@" : continue
      fields = line.split()
      if fields[2] == "*" : continue
      if fields[1] == "0": polarity = "+"
      else: polarity = "-"
      gene = fields[2]
      offset = int(fields[3])
      size = len (fields[9])
      try:
        objDic[gene].addread (polarity, offset, size)
      except KeyError:
        objDic[gene] = SmRNAwindow(gene)
        objDic[gene].addread (polarity, offset, size)
    F.close()
  else:
    print "error: invalid input format"
    sys.exit()

def z_score (table):
  value_list = [table[i] for i in sorted (table)]
  if std(value_list):
    meanlist = mean(value_list)
    stdlist = std(value_list)
    return dict ( zip ( sorted(table), [(i-meanlist)/stdlist for i in value_list] ) ) 
  else:
    return dict ( zip ( sorted(table), [0 for i in table]) )

objDic = {} # objDic is defined as a global variable
  
def __main__():
  if len(sys.argv) < 7:
    print "error: not enough parameters provided"
    sys.exit()
  load_input (sys.argv[1]) # feeds the global variable objDic (a dictionary of SmRNAwindow instances, keys=genes)
  minsize = int(sys.argv[2])
  maxsize = int(sys.argv[3])
  minscope = int(sys.argv[4])
  maxscope = int(sys.argv[5]) + 1
  general_pairs_table = dict ([(i,0) for i in range(minscope,maxscope)])
  general_prob_table = dict ([(i,0) for i in range(minscope,maxscope)])
  readcount_dic = {} # for normalized summing of local_percent_table(s)
  Total_read_in_objDic = 0
  for item in objDic:
    readcount_dic[item] = objDic[item].readcount(minsize, maxsize)
    Total_read_in_objDic += readcount_dic[item]
  OUT = open (sys.argv[-1], "w")
  for x in (objDic):
    local_pairs_table = objDic[x].count_pairs ( minsize, maxsize, range(minscope,maxscope) )
    local_prob_table = objDic[x].overlap_probability ( minsize, maxsize, range(minscope,maxscope) )
    try:
      for overlap in local_pairs_table.keys():
        general_pairs_table[overlap] = general_pairs_table.get(overlap, 0) + local_pairs_table[overlap]
    except:
      pass
    try:
      for overlap in local_prob_table.keys():
        general_prob_table[overlap] = general_prob_table.get(overlap, 0) + (1./Total_read_in_objDic*readcount_dic[x]*local_prob_table[overlap])
    except:
      pass
  z_table = z_score(general_pairs_table)
  print >> OUT, "overlap\tpairs\tz-score\toverlap_prob"
  for overlap in sorted(general_prob_table):
    print >> OUT, "%i\t%i\t%f\t%f" % (overlap, general_pairs_table[overlap], z_table[overlap], general_prob_table[overlap])
  OUT.close()

if __name__ == "__main__" : __main__()
