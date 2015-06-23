#!/usr/bin/python
# script for histogramming size of forward and reverse reads sizes in a bowtie output
# version 2
# Change in the dataframe structure for direct compatibilty with R Lattice
# Usage histogrammer_v2.py <bowtie input> <output> <Normalization Factor> <minrange> <maxrange>


import sys, re, os
from collections import defaultdict
    
def field_controle (input):
  bOUT = open (input)
  sampleline = bOUT.readline() # for test of first file line
  bOUT.close()
  samplefields = sampleline.split()
  if (samplefields[1] != "-" and samplefields[1] != "+"):
    return "error in bowtie output format"
  else:
     return "strand\tsize\tcount\titem"
  
def histogrammer (input, output, Nfactor, minrange, maxrange):
  IN = open (input)
  OUT = open (output, "w")
  Nfactor=float(Nfactor)
  print >> OUT, field_controle(input)
  hist_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int))) # hist_dic[item][forward/reverse][size 15-30nt]=count, with the enormous "autovivification" trick !!!
  for line in IN:
    fields = line.split()
    if fields[1] == "+":
      polarity = "red"
    else:
      polarity = "darkblue"
    item = fields[2]
    size = len(fields[4])
    if polarity == "red":
      hist_dict[item][polarity][size] += 1
    else:
      hist_dict[item][polarity][size] -= 1
  IN.close()
  for item in sorted(hist_dict):
    for polarity in ["red","darkblue"]:
      for size in range(minrange,maxrange+1):
        print >> OUT, "%s\t%s\t%f\t%s" % (polarity, size, hist_dict[item][polarity][size]*Nfactor, item) 
  OUT.close()


histogrammer (sys. argv[1], sys. argv[2], sys. argv[3], int(sys.argv[4]), int(sys.argv[5]) )
