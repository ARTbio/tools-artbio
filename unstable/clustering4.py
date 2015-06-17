#!/usr/bin/python
# script find clusters of small RNA reads in the genome
# version 3 - 24-12-2013 evolution to multiprocessing
# Usage clustering.py <bowtie input> <output> <bowtie index> <clustering_distance> <minimum read number per cluster to be outputed> <collapse option> <extention value> <average_cluster_size>
# <folding> <output format>


import sys, subprocess, time
from collections import defaultdict # required for some SmRNAwindow attributes (readDic)
#from numpy import mean, std # required for some SmRNAwindow methods
#from scipy import stats
from smRtools import *
import multiprocessing

def clustering (Instance):
  def clustermining (cluster, Instance): # cluster argument is a list
    if Instance.readDict[-cluster[0]]: # test whether the first position in the cluster was reverse reads
      shift = max(Instance.readDict[-cluster[0]])
      upstream_coord = cluster[0] - shift + 1
    else:
      upstream_coord = cluster[0]
    if Instance.readDict[cluster[-1]]: # test whether the last position in the cluster was forward reads
      shift = max(Instance.readDict[cluster[-1]])
      downstream_coord = cluster[-1] + shift -1
    else:
      downstream_coord = cluster[-1]
    readcount = Instance.readcount(upstream_coord=upstream_coord, downstream_coord=downstream_coord)
    mean_size, median_size, stdv_size = Instance.statsizes(upstream_coord=upstream_coord, downstream_coord=downstream_coord)
    if readcount >= minimum_reads and median_size >= min_median_size:
      location = [Instance.gene.split()[0], upstream_coord, downstream_coord]
      if output_format == "intervals":
        return "%s\t%s\t%s\t%s" % (location[0], location[1], location[2], readcount)
      cluster_size = downstream_coord - upstream_coord + 1
      if folding == "yes" and cluster_size < 151:
        foldEnergy = Instance.foldEnergy(upstream_coord=upstream_coord, downstream_coord=downstream_coord) ## be careful, test !
      else:
        foldEnergy = "."
      forwardReadcount = Instance.forwardreadcount(upstream_coord=upstream_coord, downstream_coord=downstream_coord) #
      reverseReadcount = Instance.reversereadcount(upstream_coord=upstream_coord, downstream_coord=downstream_coord) #
      density = readcount / float(cluster_size) #
      if output_format == "GFF3":
        if forwardReadcount >= reverseReadcount:
          GFFstrand = "+"
        else:
          GFFstrand = "-"
        Attributes = "ID=RC %s : FR %s : RR %s : Dens %s : Med %s : FE %s" % (readcount, forwardReadcount, reverseReadcount, density, median_size, foldEnergy)
        return "%s\tGalaxy\tRead_Cluster\t%s\t%s\t%s\t%s\t.\t%s" % (location[0], location[1], location[2], readcount, GFFstrand, Attributes)
      else:
        Forward_Barycenter, Reverse_Barycenter = Instance.barycenter(upstream_coord=upstream_coord, downstream_coord=downstream_coord)
        Zsignature = Instance.signature(24,29,24,29,range(1,27), zscore="yes", upstream_coord=upstream_coord, downstream_coord=downstream_coord)[10] #
        Hsignature = Instance.hannon_signature(24,29,24,29, range(1,27), upstream_coord=upstream_coord, downstream_coord=downstream_coord )[10] * 100
        UpiFreq = Instance.Ufreq(range(24,29), upstream_coord=upstream_coord, downstream_coord=downstream_coord)
        UsiFreq = Instance.Ufreq(range(20,22), upstream_coord=upstream_coord, downstream_coord=downstream_coord)
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (location[0], location[1], location[2], cluster_size, readcount, forwardReadcount, reverseReadcount, density, median_size, foldEnergy, Forward_Barycenter, Reverse_Barycenter, Zsignature, Hsignature, UpiFreq, UsiFreq)
    return False
  l = Instance.readDict.keys()
  l=[abs(i) for i in l]
  l=list(set(l))
  l.sort()
  upstream = 0
  cluster_list = []
  for i, element in enumerate (l[1:]):
    if abs(element-l[i]) > dist or i+2==len(l): # the 2nd part of the logical test is to capture the last cluster if it overlaps the end of the list
      cluster = l[upstream:i+1]
      upstream = i+1
      cluster_list.append(cluster)
  result_list = []
  for i in cluster_list:
    totestresult = clustermining (i, Instance)
    if totestresult: result_list.append(totestresult)
  del Instance #
  return result_list    

def logtask (results):
  global number_of_clusters
  if results:
    number_of_clusters += len(results)
    LOG.append(results)
  return
  
if __name__ == '__main__':
  start_time = time.time()
  fasta_dic = get_fasta (sys.argv[3])
  objDic = {}
  number_of_reads = 0
  F = open (sys.argv[1], "r") # F is the bowtie output taken as input
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
  
  OUT = open (sys.argv[2], "w")
  output_format=sys.argv[8]
  if output_format == "intervals":
    print >> OUT, "#chrom\tStart\tEnd\tReadCount"
  elif output_format == "GFF3":
    print >> OUT, "##gff-version 3"
  else:
    print >> OUT, "#ID\t#chrom\tStart\tEnd\tLength\tReadCount\tForwardReads\tReverseReads\tDensity\tMedian\tFoldEnergy\tForBar\tRevBar\tz-score_signature\tHannon_signature\tUfreq_in_24-28RNAs\tUfreq_in_20-21RNs"
  dist = int(sys.argv[4])
  min_median_size = int(sys.argv[6])
  minimum_reads = int(sys.argv[5])
  number_of_clusters = 0
  Instance_ID = 0
  folding=sys.argv[7]
  pool = multiprocessing.Pool(4)
  LOG = []
  instance_list = []
  for instance in objDic.keys():
    instance_list.append(objDic[instance])
  del objDic
  pool.map_async(clustering, instance_list, callback=logtask)
  pool.close()
  pool.join()
  for lines in LOG:
    for line in lines:
      print >> OUT, line
  OUT.close()
  elapsed_time = time.time() - start_time
  print "number of reads: %s\nnumber of clusters: %s\ntime: %s" % (number_of_reads, number_of_clusters, elapsed_time)
