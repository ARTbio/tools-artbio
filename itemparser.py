#!/usr/bin/python
# python parser module to analyse bowtie match whatever option nrcs or rcs has been selected
# version 2 (the previous version was only taking in charge bowtie output with +/- strand information
# Usage itemparser.py <bowtie_out> <bowtie index> <LABEL> <output file>

import sys, subprocess
from smRtools import get_fasta_headers

def SelectFilterFunction(input_file):
  F=open(input_file)
  sampleline = F.readline()
  F.close()
  samplefields=sampleline.split("\t") # bug correction for "\t" instead of any white space
  if samplefields[1] not in ["+", "-"]:
    def filter (line):
      fields = line.split("\t") # bug correction for "\t" instead of any white space
      return fields[1]
    return filter
  else:
    def filter (line):
      fields = line.split("\t") # bug correction for "\t" instead of any white space
      return fields[2]
    return filter

def produce_hitlist (bowtie_out, item_dic):
  F = open (bowtie_out, "r")
  for line in F:
    item_dic[filter(line)] += 1
  F.close()
  return item_dic

def print_hitlist (outfile, item_dic, label):
  F = open (outfile, "w")
  print >> F, "gene\t%s" % label
  for item in sorted(item_dic) :
    print >> F, "%s\t%i" % (item, item_dic[item])
  return

bowtie_out = sys.argv[1]
bowtie_index = sys.argv[2]
label = sys.argv[3]
out = sys.argv[4]
filter = SelectFilterFunction(bowtie_out)

hit_list = get_fasta_headers (bowtie_index)
for item in hit_list:
  hit_list[item] = 0

hit_list = produce_hitlist (bowtie_out, hit_list)
print_hitlist (out, hit_list, label)

