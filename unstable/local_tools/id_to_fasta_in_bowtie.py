#!/usr/bin/python

import sys

# python script to extract read id from the id column and pick up the corresponding sequence in a fasta reference read library.
# to use when 3' trimming in bowtie has degraded the information of the initial sequence read
# The script takes a bowtie standard output and create a dictionary of unique id (collapse if multimatching)
# it also takes the fasta reference read library and create another dictionary for fast matching
# it finally outputs a fasta file of reads with the original sequence (untrimmed)
# In addition, filtering can be applied to the output so that only fasta reads with the desired last three nucleotides are output.
# Usage:
# id_to_fasta_in_bowtie.py <bowtie tabular output> <fasta reference read library> <fasta with restored sequences>  <fasta or bowtie format> <filter string>

Fbowtie = open(sys.argv[1])
bowtie_dic = {}
for line in Fbowtie:
  idbowtie = line.split()[0]
  bowtie_dic[idbowtie]= line[:-1]
Fbowtie.close()

Ffasta = open(sys.argv[2])
fasta_dic = {}
for line in Ffasta:
  if line[0] == ">":
    idfasta = line[1:-1]
  else:
    fasta_dic[idfasta] = line[:-1]  
Ffasta.close()

output = open(sys.argv[3], "w")
if sys.argv[4] == "fasta": 
  try:
    sys.argv[5]
    for id in bowtie_dic:
      if fasta_dic[id][-3:]==sys.argv[5]:
        print >> output, ">%s\n%s" % (id,fasta_dic[id])
  except:
    for id in bowtie_dic:
      print >> output, ">%s\n%s" % (id,fasta_dic[id])
else:
  try:
    sys.argv[5]
    for id in bowtie_dic:
      if fasta_dic[id][-3:]==sys.argv[5]:
        fields = bowtie_dic[id].split()
        fields[3]=fasta_dic[id]
        print >> output, "\t".join(fields)
  except:
    for id in bowtie_dic:
      fields = bowtie_dic[id].split()
      fields[3]=fasta_dic[id]
      print >> output, "\t".join(fields)
output.close()
