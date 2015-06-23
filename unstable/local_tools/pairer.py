#!/usr/bin/python
# script for computing overlap signatures from a bowtie output
# version 3 17-5-2012 complete refactoring with OOP approach
# Usage pairer.py <bowtie input> <minsize query> <maxsize query> <minsize target> <maxsize target> <output> <output1> <<output2> <<output3> <<output4> <<output5> <<output6> <<bowtie index>

import sys, subprocess
from collections import defaultdict
from smRtools import get_fasta, antipara, RNAtranslate, SmRNAwindow

fasta_dic = get_fasta (sys.argv[13])
objDic = {}
F = open (sys.argv[1], "r") # F is the bowtie output taken as input

for line in F:
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

OUT = open (sys.argv[6], "w")


for x in objDic:
  sequence_list= objDic[x].pairer( 10, sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
  if sequence_list:
    for seq in sequence_list:
      print >> OUT, seq

F24nt = open (sys.argv[7], "w")
R24nt = open (sys.argv[8], "w")
F25nt = open (sys.argv[9], "w")
R25nt = open (sys.argv[10], "w")
F26nt = open (sys.argv[11], "w")
R26nt = open (sys.argv[12], "w")

n=0
for seq in sequence_list:
  n+=1
  if seq[0]=="+" and len(seq[1:]) == 24: print >> F24nt, ">%s\n%s" % (n,seq[1:])
  if seq[0]=="-" and len(seq[1:]) == 24: print >> R24nt, ">%s\n%s" % (n,seq[1:])
  if seq[0]=="+" and len(seq[1:]) == 25: print >> F25nt, ">%s\n%s" % (n,seq[1:])
  if seq[0]=="-" and len(seq[1:]) == 25: print >> R25nt, ">%s\n%s" % (n,seq[1:])
  if seq[0]=="+" and len(seq[1:]) == 26: print >> F26nt, ">%s\n%s" % (n,seq[1:])
  if seq[0]=="-" and len(seq[1:]) == 26: print >> R26nt, ">%s\n%s" % (n,seq[1:]) 

OUT.close()
F24nt.close()
R24nt.close()
F25nt.close()
R25nt.close()
F26nt.close()
R26nt.close()

base_usage = defaultdict (dict)
for i in range (23, 29):
  base_usage[i]["1U"] = 0
  base_usage[i]["10A"] = 0

for seq in sequence_list :
  if seq[1] == "U": base_usage[len(seq)-1]["1U"] += 1
  if seq[10] == "A": base_usage[len(seq)-1]["10A"] += 1

print "\t23\t24\t25\t26\t27\t28"
print "1U\t%s\t%s\t%s\t%s\t%s\t%s" % tuple (base_usage[i]["1U"] for i in range (23,29) )
print "10A\t%s\t%s\t%s\t%s\t%s\t%s" % tuple (base_usage[i]["10A"] for i in range (23,29) ) 
