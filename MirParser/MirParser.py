#!/usr/bin/python
# python parser module for pre-mir and mature miRNAs, guided by mirbase.org GFF3
# version 0.0.9 (1-6-2014)
# Usage MirParser.py  <1:index source> <2:extraction directive> <3:output pre-mir> <4: output mature miRs> <5:mirbase GFF3> <6:7:8 filePath:FileExt:FileLabel> <.. ad  lib>

import sys
from smRtools import *

IndexSource = sys.argv[1]
ExtractionDirective = sys.argv[2]
if ExtractionDirective == "--do_not_extract_index":
  genomeRefFormat = "fastaSource"
elif  ExtractionDirective == "--extract_index":
  genomeRefFormat = "bowtieIndex"
OutputPre_mirs = sys.argv[3]
OutputMature_Mirs = sys.argv[4]
GFF3_file = sys.argv[5]
Triplets = [sys.argv[6:][i:i+3] for i in xrange(0, len(sys.argv[6:]), 3)]
MasterListOfGenomes = {}

for [filePath, FileExt, FileLabel] in Triplets:
  MasterListOfGenomes[FileLabel] = HandleSmRNAwindows (filePath, FileExt, IndexSource, genomeRefFormat) 

header = ["gene"]
for [filePath, FileExt, FileLabel] in Triplets:
  header.append(FileLabel)

hit_table = ["\t".join(header)] # table header: gene, sample1, sample2, sample3, etc. separated by tabulation

## read GFF3 to subinstantiate
gff3 = open (GFF3_file, "r")
for line in gff3:
  if line[0] == "#": continue
  gff_fields = line[:-1].split("\t")
  chrom = gff_fields[0]
  gff_name = gff_fields[-1].split("Name=")[-1].split(";")[0] # to isolate the GFF Name
  item_upstream_coordinate = int(gff_fields[3])
  item_downstream_coordinate = int(gff_fields[4])
  item_polarity = gff_fields[6]
  item_line = [gff_name]
  for sample in header[1:]:
    if item_polarity == "+":
      count = MasterListOfGenomes[sample].instanceDict[chrom].forwardreadcount(upstream_coord=item_upstream_coordinate, downstream_coord=item_downstream_coordinate)
    else:
      count = MasterListOfGenomes[sample].instanceDict[chrom].reversereadcount(upstream_coord=item_upstream_coordinate, downstream_coord=item_downstream_coordinate)
    item_line.append(str(count))
  hit_table.append("\t".join(item_line) )
gff3.close()

Fpremirs = open (OutputPre_mirs, "w")
print >> Fpremirs, hit_table[0]
finalPreList = [ i for i in sorted(hit_table[1:]) if ("5p" not in i) and  ("3p" not in i)] 
print >> Fpremirs, "\n".join(finalPreList )
Fpremirs.close()

Fmaturemires = open (OutputMature_Mirs, "w")
print >> Fmaturemires, hit_table[0]
finalMatureList = [ i for i in sorted(hit_table[1:]) if ("5p" in i) or ("3p" in i)]
print >> Fmaturemires, "\n".join(finalMatureList )
Fmaturemires.close()



