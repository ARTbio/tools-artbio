#!/usr/bin/python
# script to extract cluster sequences from brennecke et al table
# version 1 - 22-03-2013
# Usage piRNA_cluster_extracter.py  <Dmfasta_genome> <coordinate reference table>

import sys

fasta_dic = {}

F = open (sys.argv[1], "r") # F the Dm fasta genome
for line in F:
  if line[0] == ">":
    try:
      fasta_dic[current_item] = "".join(stringlist) # to dump the sequence of the previous item - try because of the keyerror of the first item
    except: pass
    current_item = line[1:].rstrip().split()[0] #take the first word before space because bowtie splits headers !
    fasta_dic[current_item] = ""
    stringlist=[]
  else:
    stringlist.append(line.rstrip() )
fasta_dic[current_item] = "".join(stringlist) # for the last item
F.close()

F = open (sys.argv[2], "r") # F the cluster reference coordinates
ref_lines = F.readlines()
for line in ref_lines[1:]:
  fields = line.split()
  print ">%s_%s" % (fields[0], fields[3])
  sequence=fasta_dic[fields[0]][int(fields[1]):int(fields[2])+1]
  for i in range(0, len(sequence), 60):
    print "%s" % (sequence[i:i+60])

