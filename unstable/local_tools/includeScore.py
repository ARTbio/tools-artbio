#!/usr/bin/env python
#By, drosofff@gmail.com
# command: includeScore.py <tabulated_score_file> <GFF> <score_enriched_GFF> id_column score_column filter_limit

import sys

def get_foldchanges (file, id_column=1, score_column=5, skipline=200, limit=1E-5):
  scoredic = {} # a dic of scores with gene_id as keys
  id_column -= 1
  score_column -= 1
  linecount = 0
  F = open (file, "r")
  for line in F:
    if line[:4] == "gene" : continue
    fields = line[:-1].split("\t")
    ID = fields[id_column]
    if fields[-1]=="NA" or float(fields[-1]) > limit: continue
    try: score = float(fields[score_column])
    except: score = 0
    scoredic[ID] = score
    linecount += 1
  F.close()
  return scoredic

scoredic = get_foldchanges (sys.argv[1], limit=float(sys.argv[-1]) )


F_in = open (sys.argv[2], "r")
F_out = open (sys.argv[3], "w")


for line in F_in:
  if line[0] == "#":
    print >> F_out, line[:-1]
  else:
    GFF_fields = line[:-1].split("\t")
    GFF_attribute = GFF_fields[8]
    GFF_id = GFF_attribute.split(";")[0].split("ID=")[-1]
    print GFF_id
    if GFF_id in scoredic.keys():
      GFF_fields[5] = str(scoredic[GFF_id])
      print >> F_out, "\t".join(GFF_fields)
F_in.close()
F_out.close()

    
