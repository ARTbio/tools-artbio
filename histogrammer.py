#!/usr/bin/python
# script for histogramming size of forward and reverse reads sizes in a bowtie output
# version 1
# Usage histogrammer.py <bowtie input> <output> <Normalization Factor>


import sys, re, os
    
def field_controle (input):
  bOUT = open (input)
  sampleline = bOUT.readline() # for test of first file line
  bOUT.close()
  samplefields = sampleline.split()
  if (samplefields[1] != "-" and samplefields[1] != "+"):
    return "error in bowtie output format"
  else:
     return "size\tall\tforward\treverse"
  
def histogrammer (input, output, Nfactor):
  IN = open (input)
  OUT = open (output, "w")
  Nfactor=float(Nfactor)
  print >> OUT, field_controle(input)
  F_histogram = dict((i,0) for i in range(0,500))
  R_histogram = dict((i,0) for i in range(0,500))
  for line in IN:
    fields = line.split()
    size = len(fields[4])
    if fields[1] == "+":
      F_histogram[size] += 1
    else:
      R_histogram[size] += 1
  IN.close()
  for k in sorted(F_histogram):
#    if (k<15 or k>30): continue
    print >> OUT, "%i\t%.1f\t%.1f\t-%.1f" % (k, (F_histogram[k]+R_histogram[k])*Nfactor, F_histogram[k]*Nfactor, R_histogram[k]*Nfactor)
  OUT.close()
  return (F_histogram, R_histogram)


histogrammer (sys. argv[1], sys. argv[2], sys. argv[3])
