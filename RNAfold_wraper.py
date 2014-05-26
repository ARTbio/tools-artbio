#!/usr/bin/python
import sys
from smRtools import *

F=open (sys.argv[1])
crudelist=F.read()
lines=crudelist.split("\n")
lines=lines[1:-1]
sequences=[]
for line in lines:
  sequence=line.split()[-1]
  sequences.append(sequence)
print sequences
energies= RNAfold(sequences)
for i in energies:
  print "%s\t%s" % (i, energies[i])
