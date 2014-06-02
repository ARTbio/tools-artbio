#!/usr/bin/python
# python parser module to analyse sRbowtie alignments
# version 0.9
# Usage sRbowtieParser.py  <1:index source> <2:extraction directive> <3:outputL> <4:polarity> <5:6:7 filePath:FileExt:FileLabel> <.. ad  lib>

import sys
from smRtools import *

IndexSource = sys.argv[1]
ExtractionDirective = sys.argv[2]
if ExtractionDirective == "--do_not_extract_index":
  genomeRefFormat = "fastaSource"
elif  ExtractionDirective == "--extract_index":
  genomeRefFormat = "bowtieIndex"
Output = sys.argv[3]
Polarity = sys.argv[4] # maybe "both", "sense", "antisense"
Triplets = [sys.argv[5:][i:i+3] for i in xrange(0, len(sys.argv[5:]), 3)]
MasterListOfGenomes = {}

for [filePath, FileExt, FileLabel] in Triplets:
  MasterListOfGenomes[FileLabel] = HandleSmRNAwindows (filePath, FileExt, IndexSource, genomeRefFormat) 

header = ["gene"]
for [filePath, FileExt, FileLabel] in Triplets:
  header.append(FileLabel)

F = open (sys.argv[3], "w")
print >> F, "\t".join(header)
for item in sorted (MasterListOfGenomes[header[1]].instanceDict.keys() ):
  line=[item]
  for sample in header[1:]:
    if Polarity == "both":
      count = str (MasterListOfGenomes[sample].instanceDict[item].readcount())
    elif Polarity == "sense":
      count = str (MasterListOfGenomes[sample].instanceDict[item].forwardreadcount())
    elif Polarity == "antisense":
      count = str (MasterListOfGenomes[sample].instanceDict[item].reversereadcount())
    line.append(count)
  print >> F,  "\t".join(line )
F.close()
