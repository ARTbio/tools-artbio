#!/usr/bin/python
# script for generating 4 fasta files 21nt and 25nt, forward and reverse, before weblogo analysis
# version 3-2-2012
# Usage fasta_sorter.py <bowtie input> <output1> <output2> <output3> <output4>

import sys, re, os

def antipara (sequence):
  antidict = {"A":"T", "T":"A", "G":"C", "C":"G"}
  revseq = sequence[::-1]
  return "".join([antidict[i] for i in revseq])

def RNAtranslate (sequence):
  return "".join([i if i in "AGC" else "U" for i in sequence])

def dispatch (bowtie_input, f21, r21, f25, r25):
  IN = open (bowtie_input)
  F21= open (f21, "w")
  R21= open (r21, "w")
  F25= open (f25, "w")
  R25= open (r25, "w")
  for line in IN:
    fields = line.split()
    read_header = fields[0]
    read_polarity = fields[1]
    read_sequence = fields[4]
    if "N" in read_sequence: continue
    read_size = len(read_sequence)
    if read_polarity == "+" and read_size == 25:
      seq = RNAtranslate (read_sequence)
      print >> F25, ">%s\n%s" % (read_header, seq)
    elif read_polarity == "-" and read_size == 25:
       seq = RNAtranslate (antipara(read_sequence))
       print >> R25, ">%s\n%s" % (read_header, seq)    
    elif read_polarity == "+" and read_size == 21:
      seq = RNAtranslate (read_sequence)
      print >> F21, ">%s\n%s" % (read_header, seq)
    elif read_polarity == "-" and read_size == 21:
       seq = RNAtranslate (antipara(read_sequence))
       print >> R21, ">%s\n%s" % (read_header, seq)
  IN.close()
  F21.close()
  R21.close()
  F25.close()
  R25.close()
  return  

dispatch (sys. argv[1], sys. argv[2], sys. argv[3], sys. argv[4], sys. argv[5])



