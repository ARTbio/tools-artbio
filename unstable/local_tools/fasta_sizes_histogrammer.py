#!/usr/bin/python
# script for histogramming size of sequences in a number of fasta files
# version 1 27-9-2013
# Usage fasta_sizes_histogrammer.py <file1> <label1> <file2> <label2> ... <output>


import sys, re, os
    
def histogrammer (file, label):
  IN = open (file)
  histogram = dict((i,0) for i in range(0,102))
  for line in IN:
    if line[0] == ">": continue
    else:
      histogram[len(line[:-1])] += 1
  IN.close()
  result_table = []
  for k in sorted(histogram):
    result_table.append ( "%i\t%i\t%s" % (k, histogram[k], label) )
  return (result_table)

F = open(sys.argv[-1], "w")
print >> F, "size\treads\tsample"
for file, label in zip (sys.argv[1:-1:2], sys.argv[2:-1:2]) :
  filedata = histogrammer (file, label)
  for line in filedata:
    print >> F, line
F.close()  
