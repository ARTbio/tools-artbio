#!/usr/bin/python
# python parser module to analyse sRbowtie alignments
# version 1.0.2 - argparse implementation
# Usage sRbowtieParser.py  <1:index source> <2:extraction directive> <3:outputL> <4:polarity> <5:6:7 filePath:FileExt:FileLabel> <.. ad  lib>

import sys, argparse
from smRtools import *

def masterListGenerator(data_source):
  for filePath, FileExt, FileLabel in data_source:
    yield HandleSmRNAwindows (filePath, FileExt, IndexSource, genomeRefFormat)

def Parser():
  the_parser = argparse.ArgumentParser()
  the_parser.add_argument('--IndexSource', action="store", type=str, help="Path to the index source")
  the_parser.add_argument('--ExtractDirective', action="store", type=str, choices=["fastaSource", "bowtieIndex"], help="Extract info from fasta or bowtie index")
  the_parser.add_argument('--output', action="store", type=str, help="path to the output")
  the_parser.add_argument('--polarity', choices=["forward", "reverse", "both"], help="forward, reverse or both forward an reverse reads are counted")
  the_parser.add_argument('--alignmentSource',nargs='+', help="paths to alignments files")
  the_parser.add_argument('--alignmentFormat',nargs='+', help="Format of the bowtie alignment (tabular, sam or bam)")
  the_parser.add_argument('--alignmentLabel',nargs='+', help="Label of the alignment")
  args = the_parser.parse_args()
  return args

args = Parser()

IndexSource = args.IndexSource
genomeRefFormat = args.ExtractDirective
Output = args.output
Polarity = args.polarity
header = ["gene"]


FileLabelList=[label for label in args.alignmentLabel]
header.extend(FileLabelList)
assert (len(FileLabelList)==len(set(FileLabelList))),"You have supplied a non-unique label. Please make sure that your input files have unique names"

data_source=zip (args.alignmentSource, args.alignmentFormat, args.alignmentLabel)
master_generator=masterListGenerator(data_source)

for i,window in enumerate(master_generator):
  window=window
  if i==0:
    gene_count_dict={gene:[str(item.readcount(polarity=Polarity))] for gene,item in window.instanceDict.items()}
  else:
    [gene_count_dict[gene].append(str(item.readcount(polarity=Polarity))) for gene,item in window.instanceDict.items()]


F = open (args.output, "w")
# print >>F, args
print >> F, "\t".join(header)

for item in sorted(gene_count_dict.keys()):
  line=[item]
  line.extend(gene_count_dict[item])
  print >> F,  "\t".join(line )
F.close()
