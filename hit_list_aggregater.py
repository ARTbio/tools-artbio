#!/usr/bin/python
# python parser module to extract miRNA hit list
# drosofff@gmail.com
# version 2 (for a flatten mir 5p and 3p hit list)
# Usage hit_list_aggregater2.py <outlist> <mir_hitlist1> <mir_hitlist2> ...

import sys, re, os, string
from collections import defaultdict
class recursivedefaultdict(defaultdict):
  def __init__(self):
    self.default_factory = type(self)
    
hash_results = recursivedefaultdict()

def extract_list (file, tag, hash_results):
  F=open(file, "r")
  Flag = False
  for line in F:
    if line.startswith(tag):
      Flag = True
      label = line[:-1].split("\t")[1]
      continue
    if Flag:
      fields = line[:-1].split("\t")
      name = fields[0]
      hash_results[name][label] = fields[1]
  F.close()
  return label
      

def format_hash (hash_results, list_of_labels, output):
  header = "gene\t" + "\t".join ( list_of_labels ) 
  F = open (output, "w")
  print >> F, header
  for mir in sorted (hash_results.keys()):
    string = mir 
    for label in list_of_labels:
      string += "\t%s" % (hash_results[mir][label])
    print >> F, string
  F.close()
  
list_of_labels = []
for file in sys.argv[2:]:
  list_of_labels.append (extract_list (file, "gene", hash_results) )
  
format_hash (hash_results, list_of_labels, sys.argv[1])

