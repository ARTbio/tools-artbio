#!/usr/bin/env python
# from http://blog.nextgenetics.net/?e=27
# usage GFF3toGTF.py <GFF3 file> <GTF converted file>

import sys

inFile = open(sys.argv[1],'r')
outFile = open(sys.argv[2],'w')

for line in inFile:
  #skip comment lines that start with the '#' character
  if line[0] != '#':
    theline = line[:-1]
    if theline.split("\t")[2] != "exon": continue
    theline = theline.replace ("gene_id", "FBgene")
    theline = theline.replace ("transcript_id", "FBtr")
    theline = theline.replace ("transcript_name", "transcript_id") 
    #print out this new GTF line
    print >> outFile, theline

inFile.close()
outFile.close()
