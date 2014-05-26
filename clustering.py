#!/usr/bin/python
# script find clusters of small RNA reads in the genome
# version 2 - 22-11-2013 evolution to use a subinstance extraction approach
# Usage clustering.py <bowtie input> <output> <bowtie index> <clustering_distance> <minimum read number per cluster to be outputed> <collapse option> <extention value> <average_cluster_size>
# <folding> <output format>


import sys, subprocess
from collections import defaultdict # required for some SmRNAwindow attributes (readDic)
#from numpy import mean, std # required for some SmRNAwindow methods
from scipy import stats
from smRtools import *

def clustermining (cluster, Instance_ID): # cluster argument is a list
  if objDic[object].readDict[-cluster[0]]: # test whether the first position in the cluster was reverse reads
    shift = max(objDic[object].readDict[-cluster[0]])
    upstream_coord = cluster[0] - shift + 1
  else:
    shift = max(objDic[object].readDict[cluster[0]])
    upstream_coord = cluster[0]
  if objDic[object].readDict[cluster[-1]]: # test whether the last position in the cluster was forward reads
    shift = max(objDic[object].readDict[cluster[-1]])
    downstream_coord = cluster[-1] + shift -1
  else:
    shift = max(objDic[object].readDict[-cluster[-1]])
    downstream_coord = cluster[-1]
  subinstanceDic = extractsubinstance (upstream_coord, downstream_coord, objDic[object])
  if subinstanceDic.readcount() >= minimum_reads and subinstanceDic.statsizes()[1] > min_median_size:
    global number_of_clusters
    number_of_clusters += 1
    location = subinstanceDic.gene.split()
    cluster_size = subinstanceDic.size
    if folding == "yes" and cluster_size < 151:
      foldEnergy = subinstanceDic.foldEnergy()
    else:
      foldEnergy = "."
    sequence = subinstanceDic.sequence
    readcount = subinstanceDic.readcount()
    forwardReadcount = subinstanceDic.forwardreadcount()
    reverseReadcount = subinstanceDic.reversereadcount()
    density = readcount / float(cluster_size)
    stdv = subinstanceDic.statsizes()[2]
    median = subinstanceDic.statsizes()[1]
    if output_format == "intervals":
      print >> OUT, "%s\t%s\t%s\t%s" % (location[0], location[1], location[2], readcount)
    elif output_format == "GFF3":
      if forwardReadcount >= reverseReadcount:
        GFFstrand = "+"
      else:
        GFFstrand = "-"
      Attributes = "ID=RC %s : FR %s : RR %s : Dens %s : Med %s : FE %s" % (readcount, forwardReadcount, reverseReadcount, density, median, foldEnergy)
      print >> OUT, "%s\tGalaxy\tRead_Cluster\t%s\t%s\t%s\t%s\t.\t%s" % (location[0], location[1], location[2], readcount, GFFstrand, Attributes)
    else:
      Forward_Barycenter = subinstanceDic.barycenter()[0]
      Reverse_Barycenter = subinstanceDic.barycenter()[1]
      Zsignature = subinstanceDic.signature(24,29,24,29,range(1,27), zscore="yes" )[10]
      Hsignature = subinstanceDic.hannon_signature(24,29,24,29, range(1,27) )[10] * 100
      UpiFreq = subinstanceDic.Ufreq(range(24,29))
      UsiFreq = subinstanceDic.Ufreq(range(20,22))
      print >> OUT, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (Instance_ID, location[0], location[1], location[2], cluster_size, readcount, forwardReadcount, reverseReadcount, density, median, foldEnergy, Forward_Barycenter, Reverse_Barycenter, Zsignature, Hsignature, UpiFreq, UsiFreq)
  del subinstanceDic
  return

OUT = open (sys.argv[2], "w")
dist = int(sys.argv[4])
min_median_size = int(sys.argv[6])
minimum_reads = int(sys.argv[5])
fasta_dic = get_fasta (sys.argv[3])
folding=sys.argv[7]
output_format=sys.argv[8]
objDic = {}
F = open (sys.argv[1], "r") # F is the bowtie output taken as input
number_of_reads = 0
number_of_clusters = 0

for line in F:
  number_of_reads += 1
  fields = line.split()
  polarity = fields[1]
  gene = fields[2]
  offset = int(fields[3])
  size = len (fields[4])
  try:
    objDic[gene].addread (polarity, offset, size)
  except KeyError:
    objDic[gene] = SmRNAwindow(gene, fasta_dic[gene])
    objDic[gene].addread (polarity, offset, size)
F.close()
Instance_ID = 0
if output_format == "intervals":
  print >> OUT, "#chrom\tStart\tEnd\tReadCount"
elif output_format == "GFF3":
  print >> OUT, "##gff-version 3"
else:
  print >> OUT, "#ID\t#chrom\tStart\tEnd\tLength\tReadCount\tForwardReads\tReverseReads\tDensity\tMedian\tFoldEnergy\tForBar\tRevBar\tz-score_signature\tHannon_signature\tUfreq_in_24-28RNAs\tUfreq_in_20-21RNs"

stupidlist=objDic.keys()

for object in stupidlist:
  l = objDic[object].readDict.keys()
  l=[abs(i) for i in l]
  l=list(set(l))
  l.sort()
  upstream = 0
  for i, element in enumerate (l[1:]):
    if abs(element-l[i]) > dist:
      cluster = l[upstream:i+1]
      upstream = i+1
      Instance_ID += 1
      clustermining (cluster, Instance_ID)
  Instance_ID += 1
  clustermining (l[upstream:], Instance_ID) # dernier cluster ? to test
  del objDic[object]
OUT.close()

print "number of reads: %s\nnumber of clusters: %s" % (number_of_reads, number_of_clusters)
