#!/usr/bin/python
# python parser module for for readmaps and size distributions, guided by GFF3
# version 0.9.1 (1-6-2014)
# Usage readmap.py  <1:index source> <2:extraction directive> <3:output pre-mir> <4: output mature miRs> <5:mirbase GFF3>
#                     <6:pathToLatticeDataframe or "dummy_dataframe_path"> <7:Rcode or "dummy_plotCode"> <8:latticePDF or "dummy_latticePDF">
#                     <9:10:11 filePath:FileExt:FileLabel> <.. ad  lib>

import sys, subprocess, argparse
from smRtools import *
from collections import OrderedDict, defaultdict
import os

def Parser():
  the_parser = argparse.ArgumentParser()
  the_parser.add_argument('--output_readmap', action="store", type=str, help="readmap dataframe")
  the_parser.add_argument('--output_size_distribution', action="store", type=str, help="size distribution dataframe")
  the_parser.add_argument('--reference_fasta', action="store", type=str, help="output file")
  the_parser.add_argument('--reference_bowtie_index',action='store', help="paths to indexed or fasta references")
  the_parser.add_argument('--input',nargs='+', help="paths to multiple input files")
  the_parser.add_argument('--ext',nargs='+', help="input file type")
  the_parser.add_argument('--label',nargs='+', help="labels of multiple input files")
  the_parser.add_argument('--normalization_factor',nargs='+', type=float, help="Normalization factor for input file")
  the_parser.add_argument('--gff', type=str, help="GFF containing regions of interest")
  the_parser.add_argument('--minquery', type=int, help="Minimum readsize")
  the_parser.add_argument('--maxquery', type=int, help="Maximum readsize")
  the_parser.add_argument('--rcode', type=str, help="R script")
  args = the_parser.parse_args()
  return args

args=Parser()
if args.reference_fasta:
  genomeRefFormat = "fastaSource"
  genomeRefFile = args.reference_fasta  
if args.reference_bowtie_index:
  genomeRefFormat = "bowtieIndex"
  genomeRefFile = args.reference_bowtie_index  
readmap_file=args.output_readmap
size_distribution_file=args.output_size_distribution
minquery=args.minquery
maxquery=args.maxquery
Rcode = args.rcode
filePath=args.input
fileExt=args.ext
fileLabel=args.label
normalization_factor=args.normalization_factor

MasterListOfGenomes = OrderedDict()

def process_samples(filePath):
  for i, filePath in enumerate(filePath):
    norm=normalization_factor[i]
    print fileLabel[i]
    MasterListOfGenomes[fileLabel[i]] = HandleSmRNAwindows (alignmentFile=filePath, alignmentFileFormat=fileExt[i], genomeRefFile=genomeRefFile, genomeRefFormat=genomeRefFormat,\
                        biosample=fileLabel[i], size_inf=minquery, size_sup=maxquery, norm=norm)
  return MasterListOfGenomes

def dataframe_sanityzer (listofdatalines):
  Dict = defaultdict(float) 
  for line in listofdatalines:
    fields= line.split("\t")
    Dict[fields[0]] += float (fields[2])
  filtered_list = []
  for line in listofdatalines:
    fields= line.split("\t")
    if Dict[fields[0]] != 0:
      filtered_list.append(line) 
  return filtered_list


def listify_plottable_item(item):
  """
  plottable is a list of strings:
  'FBti0020401\t78\t-1.0\tR'
  split on tab and return gene, coordinate, count and orientation
  """
  gene, coordinate, count, orientation = item.split("\t")
  return gene, coordinate, count, orientation

def lookup_gene_length(gene, readDict):
  return readDict[readDict.keys()[0]].instanceDict[gene].size

def handle_start_stop_coordinates(plottable, readDict):
  """
  To ensure that the plot area always includes the correct start and end coordinates,
  we add an entry at start [coordinate 0] and end [last coordinate] of count 0, if these do not exist.
  """
  first_line = plottable[0]
  last_line = plottable[-1]
  gene, coordinate, count, orientation = listify_plottable_item(first_line)
  if not coordinate == "0":
    new_line = "\t".join([gene, "0", "0", "F"])
    plottable = [new_line] + plottable
  gene_length = str(lookup_gene_length(gene, readDict))
  gene, coordinate, count, orientation = listify_plottable_item(last_line)
  if not coordinate == gene_length:
    last_line = "\t".join([gene, gene_length, "0", "F"])
    plottable = plottable + [last_line]
  return plottable

def write_readplot_dataframe(readDict, readmap_file):
  listoflines = []
  with open(readmap_file, 'w') as readmap:
    print >>readmap, "gene\tcoord\tcount\tpolarity\tsample"
    for sample in readDict.keys():
      if args.gff:
        dict=readDict[sample]
      else:
        dict=readDict[sample].instanceDict
      for gene in dict.keys():
        plottable = dict[gene].readplot()
        plottable = handle_start_stop_coordinates(plottable, readDict)
        for line in plottable:
          #print >>readmap, "%s\t%s" % (line, sample)
          listoflines.append ("%s\t%s" % (line, sample))
    listoflines = dataframe_sanityzer(listoflines)
    for line in listoflines:
      print >>readmap, line

def write_size_distribution_dataframe(readDict, size_distribution_file):
  listoflines = []
  with open(size_distribution_file, 'w') as size_distrib:
    print >>size_distrib, "gene\tsize\tcount\tpolarity\tsample" # test before was "gene\tpolarity\tsize\tcount\tsample"
    for sample in readDict.keys():
      if args.gff:
        dict=readDict[sample]
      else:
        dict=readDict[sample].instanceDict
      for gene in dict.keys():
        histogram = dict[gene].size_histogram(minquery=args.minquery, maxquery=args.maxquery)
        for polarity in histogram.keys():
          if polarity=='both':
            continue
          #for size in xrange(args.minquery, args.maxquery):
          #  if not size in histogram[polarity].keys():
          #    histogram[size]=0
          for size, count in histogram[polarity].iteritems():
            #print >>size_distrib, "%s\t%s\t%s\t%s\t%s" % (gene, size, count, polarity, sample) # test, changed the order accordingly
            listoflines.append ("%s\t%s\t%s\t%s\t%s" % (gene, size, count, polarity, sample) )
    listoflines = dataframe_sanityzer(listoflines)
    for line in listoflines:
      print >>size_distrib, line  

def gff_item_subinstances(readDict, gff3):
  GFFinstanceDict=OrderedDict()
  for sample in readDict.keys():
    GFFinstanceDict[sample]={} # to implement the 2nd level of directionary in an OrderedDict Class object (would not be required with defaultdict Class)
  with open(gff3) as gff:
    for line in gff:
      if line[0] == "#": continue
      gff_fields = line[:-1].split("\t")
      chrom = gff_fields[0]
      gff_name = gff_fields[-1].split("Name=")[-1].split(";")[0] # to isolate the GFF Name
      item_upstream_coordinate = int(gff_fields[3])
      item_downstream_coordinate = int(gff_fields[4])
      item_polarity = gff_fields[6]
      for sample in readDict.keys():
## this is not required anymore but test
#        if not GFFinstanceDict.has_key(sample):
#          GFFinstanceDict[sample]={}
####
        subinstance=extractsubinstance(item_upstream_coordinate, item_downstream_coordinate, readDict[sample].instanceDict[chrom])
        if item_polarity == '-':
          subinstance.readDict={key*-1:value for key, value in subinstance.readDict.iteritems()}
        subinstance.gene=gff_name
        GFFinstanceDict[sample][gff_name]=subinstance
  return GFFinstanceDict

MasterListOfGenomes=process_samples(filePath)

if args.gff:
  MasterListOfGenomes=gff_item_subinstances(MasterListOfGenomes, args.gff)

write_readplot_dataframe(MasterListOfGenomes, readmap_file)
write_size_distribution_dataframe(MasterListOfGenomes, size_distribution_file)

R_command="Rscript "+ Rcode
process = subprocess.Popen(R_command.split())
process.wait()
	
