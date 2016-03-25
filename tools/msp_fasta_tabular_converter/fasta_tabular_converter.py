#!/usr/bin/python
#
import sys
import string
import argparse
from collections import defaultdict

def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="input file")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output converted file")
    the_parser.add_argument(
        '--type', action="store", type=str, help="type of convertion")
    args = the_parser.parse_args()
    return args

def readfasta_writetabular(fasta, tabular, mode="oneline"):
    F = open(fasta, "r")
    for line in F:
        if line[0] == ">":
            try:
                seqdic["".join(stringlist)] += 1 # to dump the sequence of the previous item - try because of first missing stringlist variable
            except: pass
            stringlist=[]
        else:
            stringlist.append(line[:-1])
    try:
        seqdic["".join(stringlist)] +=  1 # for the last sequence
    except: pass # in case file to convert is empty
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

def main(input, output, type):
    if type == "fasta2tabular":
        readfasta_writetabular(input, output)
    elif type == "tabular2fasta":
        readtabular_writefasta(input, output)
    elif type == "tabular2fastaweight":
        readtabular_writefastaweighted (input, output)
    elif type == "fastaweight2fastaweight":
        readfastaeighted_writefastaweighted(input, output)
    elif type == "fastaweight2fasta":
        readfastaeighted_writefasta(input, output)

if __name__ == "__main__":
    seqdic = defaultdict(int)
    args = Parser()
    main (args.input, args.output, args.type)
