#!/usr/bin/python
# script extract a particular genome region
# version 3 - 4-12-2012 to simplify the previous 27-12-2011 version with "usine a gaz" class
# Usage bowtie_window_extraction.py <bowtie input> <geneID> <Upstream_coordinate> <Downstream_coordinate> <output> <minsize> <maxsize> <output_format>


import sys, subprocess
from smRtools import antipara


geneID = sys.argv[2]
Upstream_coordinate = int(sys.argv[3])
Downstream_coordinate = int(sys.argv[4])
minsize = int(sys.argv[6])
maxsize = int(sys.argv[7])
F = open (sys.argv[1], "r") # F is the bowtie output taken as input
OUT = open (sys.argv[5], "w")
output_format = sys.argv[8]

def selectorformat(format, file_handler):
  if format == "bowtie":
    print "bowtie format"
    def printformat(header, polarity, geneID, coordinate, sequence, file_handler):
      print >> file_handler, "%s\t%s\t%s\t%s\t%s" % (header, polarity, geneID, coordinate, sequence)
  elif format == "fasta":
    print "fasta format"
    def printformat(header, polarity, geneID, coordinate, sequence, file_handler):
      if polarity == "-": sequence = antipara(sequence)
      print >> file_handler, ">%s\n%s" % (header, sequence)
  return printformat

#printing format selection
theprint = selectorformat(output_format, OUT)

def selectorfilter(option):
  if option == "by_location":
    print "choosing by_location"
    def filter(header, polarity, running_geneID, running_coor, sequence, file_handler):
      if geneID != running_geneID:
        return
      if polarity == "-":
        original_coor = running_coor
        running_coor = running_coor + len(sequence) -1
      else: original_coor = running_coor
      if Upstream_coordinate <=running_coor <=  Downstream_coordinate :
        theprint (header, polarity, running_geneID, original_coor, sequence, file_handler)

  if option == "by_size":
    print "choosing by_size"
    def filter(header, polarity, running_geneID, running_coor, sequence, file_handler):
      running_size = len(sequence)
      if minsize <= running_size <= maxsize:
        theprint (header, polarity, running_geneID, running_coor, sequence, file_handler)

  if option == "by_location&size":
    print "choosing by_location&size" 
    def filter(header, polarity, running_geneID, running_coor, sequence, file_handler):
      if geneID != running_geneID:
        return
      if polarity == "-":
        original_coor = running_coor
        running_coor = running_coor + len(sequence) -1
      else:
        original_coor = running_coor
      if not(Upstream_coordinate <= running_coor <= Downstream_coordinate):
        return
      running_size = len (sequence)
      if not(minsize <= running_size <= maxsize):
        return
      theprint (header, polarity, running_geneID, original_coor, sequence, file_handler)
  return filter


### filter selection
if geneID!="item" and maxsize:
  thefilter=selectorfilter("by_location&size")
elif geneID!="item":
  thefilter=selectorfilter("by_location")
else:
  thefilter=selectorfilter("by_size")

### filtering
for line in F:
  fields = line.split()
  header = fields[0]
  polarity = fields[1]
  running_geneID = fields[2]
  coordinate = int(fields[3])+1 # to shift to 1-based coordinates of genome browser and humans ! (bowtie is 0-based)
  sequence = fields[4]
  thefilter(header, polarity, running_geneID, coordinate, sequence, OUT) 

F.close()
OUT.close()
