#!/usr/bin/python
#
import sys
from collections import defaultdict

def readfasta_writetabular(fasta, tabular):
  F = open(fasta, "r")
  for line in F:
    if line[0] == ">": continue
    else:
      seqdic[line[:-1]] += 1
  F.close()
  F = open(tabular, "w")
  for seq in sorted(seqdic, key=seqdic.get, reverse=True):
    print >> F, "%s\t%s" % (seq, seqdic[seq])
  F.close()
    
        
def readtabular_writefasta(tabular, fasta):
  F = open(tabular, "r")
  Fw = open(fasta, "w")
  counter = 0
  for line in F:
    fields = line.split()
    for i in range(int(fields[1])):
      counter += 1
      print >> Fw, ">%s\n%s" % (counter, fields[0])
  F.close()
  Fw.close()

def readtabular_writefastaweighted (tabular, fasta):
  F = open(tabular, "r")
  Fw = open(fasta, "w")
  counter = 0
  for line in F:
    counter += 1
    fields = line[:-1].split()
    print >> Fw, ">%s_%s\n%s" % (counter, fields[1],  fields[0])
  F.close()
  Fw.close()

def readfastaeighted_writefastaweighted(fastaweigthed_input, fastaweigthed_reparsed):
  F = open(fastaweigthed_input, "r")
  number_reads = 0
  for line in F:
    if line[0] == ">":
      weigth = int(line[1:-1].split("_")[-1])
      number_reads += weigth
    else:
      seqdic[line[:-1]] += weigth
  F.close()
  F = open(fastaweigthed_reparsed, "w")
  n=0
  for seq in sorted(seqdic, key=seqdic.get, reverse=True):
    n += 1
    print >> F, ">%s_%s\n%s" % (n, seqdic[seq], seq)
  F.close()
  print "%s reads collapsed" % number_reads

def readfastaeighted_writefasta(fastaweigthed, fasta):
  F = open(fastaweigthed, "r")
  Fw = open(fasta, "w")
  counter = 0
  for line in F:
    if line[0] == ">":
      weigth = int(line[1:-1].split("_")[-1])
    else:
      seq = line[:-1]
      for i in range (weigth):
        counter += 1
        print >> Fw, ">%s\n%s" % (counter, seq)
  F.close()
  Fw.close()


seqdic = defaultdict(int)
option = sys.argv[3]

if option == "fasta2tabular":
  readfasta_writetabular(sys.argv[1], sys.argv[2])
elif option == "tabular2fasta":
  readtabular_writefasta(sys.argv[1], sys.argv[2])
elif option == "tabular2fastaweight":
  readtabular_writefastaweighted (sys.argv[1], sys.argv[2])
elif option == "fastaweight2fastaweight":
  readfastaeighted_writefastaweighted(sys.argv[1], sys.argv[2])
elif option == "fastaweight2fasta":
  readfastaeighted_writefasta(sys.argv[1], sys.argv[2])
