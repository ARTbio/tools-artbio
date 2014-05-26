#!/usr/bin/python
# script for generating 2 fasta files of 23nt-TRUNCATED 23-28nt reads, forward and reverse, before weblogo analysis
# version 23-5-2012
# Usage trunc_pi23.py <bowtie input> <output1> <output2>

import sys, re, os

def antipara (sequence):
  antidict = {"A":"T", "T":"A", "G":"C", "C":"G"}
  revseq = sequence[::-1]
  return "".join([antidict[i] for i in revseq])

def RNAtranslate (sequence):
  return "".join([i if i in "AGC" else "U" for i in sequence])

def dispatch (bowtie_input, f23, r23):
  IN = open (bowtie_input)
  F23= open (f23, "w")
  R23= open (r23, "w")
  for line in IN:
    fields = line.split()
    read_header = fields[0]
    read_polarity = fields[1]
    read_sequence = fields[4]
    if "N" in read_sequence: continue
    read_size = len(read_sequence)
    if read_polarity == "+" and 23<read_size<28:
      seq = RNAtranslate (read_sequence)
      print >> F23, ">%s\n%s" % (read_header, seq[:23])
    elif read_polarity == "-" and 23<read_size<28:
       seq = RNAtranslate (antipara(read_sequence))
       print >> R23, ">%s\n%s" % (read_header, seq[:23])    
  IN.close()
  F23.close()
  R23.close()
  return  

dispatch (sys. argv[1], sys. argv[2], sys. argv[3])



