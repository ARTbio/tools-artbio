#!/usr/bin/python
# python parser module for pipeline for miRNA full profiling with bowtie 20/4/2011
# version 1
# Usage class_mirparser.py <bowtie_out> <bowtie miRNA index> <LABEL> <output1 file> <output2 file>
# still to work. Split count still bugged
import sys, subprocess
from collections import defaultdict

def get_fasta (index="/home/galaxy/galaxy-dist/bowtie/5.37_Dmel/5.37_Dmel"): # here one must supply the exact base path for the library in galaxy environment
  p = subprocess.Popen(args=["bowtie-inspect","-a", "0", index], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  outputlines = p.stdout.readlines()
  p.wait()
  miR_sequence_liste = {}
  for line in outputlines:
    if line[0] == ">":
      miR_name = line[1:-1].split()[0]
    else:
      miR_sequence_liste[miR_name] = line[:-1]
  return miR_sequence_liste
  
class Mirna:

  def __init__(self, name, sequence):
    self.name = name
    self.sequence = sequence
    self.matched_reads = []
    self.dicmap = {}

  def addread (self, offset, size):
    #method to add reads to the object
    self.matched_reads.append( (offset, size) )
    self.dicmap[(offset, size)] = self.dicmap.get((offset, size), 0) + 1
    return

  def mircount (self):
    #method to return the raw counts
    return len(self.matched_reads)
    
  def density (self):
    '''method to output the read coverage by position in the mir'''
    map = [0 for i in range (len(self.sequence))]
    for offset, size in self.dicmap:
      for i in range (offset, offset+size):
        map[i] += self.dicmap[(offset,size)]
    return map

  def normalized_density (self):
    map = self.density ()
    maximum = float (max (map) ) or 1
    length = float (len (map) ) or 1
    Total_NoR = self.mircount()
    output = ["mir\tcoordinate\tdensity\tNoR"]
    for i, D in enumerate (map):
      output.append("%s\t%s\t%s\t%s" % (self.name, (i+1)/length, D/maximum, Total_NoR))
    return "\n".join(output)

  def hitmap (self):
    #method to output the reads above the premir
    output = []
    output.append (self.name)
    output.append ( "%s\t%s\t%s\t%s" % (self.sequence, "offset", "size", "num reads") )
    for pos_size in sorted(self.dicmap):
      seq = self.sequence[ pos_size[0] : pos_size[0]+pos_size[1] ]
      output.append ("%s%s%s\t%s\t%s\t%s" % ("."*len(self.sequence[:pos_size[0]]), seq, "."*len(self.sequence[pos_size[0]+pos_size[1]:]), pos_size[0]+1, pos_size[1], self.dicmap[pos_size] ) ) #attention pos_size[0]+1 because 1-based offset for biologists
    return "\n".join(output)

  def splitcount (self, shift):
    #method to assign counts to 5p and 3p parts
    median = len(self.sequence)/2
    splitsite = 0
    scores = []
    for i in range(median-shift, median+shift+1):
      countsum = 0
      for pos_size in self.dicmap:
        if pos_size[0] <= i <= pos_size[0]+pos_size[1]-1: continue
        else: countsum = countsum + self.dicmap[pos_size]
      scores.append(countsum)      
    firstmax = scores.index(max(scores))
    scores.reverse()
    lastmax = scores.index(max(scores))
    scores.reverse()
    split_selected = firstmax + len(scores[firstmax:-lastmax])/2 + median - shift
    mir5p = 0
    mir3p = 0
    for pos_size in self.dicmap:
      if pos_size[0] <= split_selected <= pos_size[0]+pos_size[1]-1: continue
      elif split_selected <= pos_size[0]: mir3p = mir3p + self.dicmap[pos_size]
      else : mir5p = mir5p + self.dicmap[pos_size]
    return  "%s_5p\t%s\n%s_3p\t%s" % (self.name, mir5p, self.name, mir3p)
    
  def splitcount_2 (self, shift):
    #new method to assign counts to 5p and 3p parts, base on density map
    density_map = self.density()
    median = len(self.sequence)/2
    minimum = 0
    densitydic = dict ([(i, density) for i, density in enumerate (density_map) if median-shift<= i <= median+shift ])
    revdic = dict(map(lambda item: (item[1],item[0]),densitydic.items()))
    mindensity_offset = revdic[min(revdic.keys())]
    mir5p = 0
    mir3p = 0
    for pos_size in self.dicmap:
      if mindensity_offset in range (pos_size[0], pos_size[0]+pos_size[1]): continue
#      if pos_size[0] <= mindensity_offset <= pos_size[0]+pos_size[1]-1: continue
      if mindensity_offset <= pos_size[0]: mir3p = mir3p + self.dicmap[pos_size]
      else : mir5p = mir5p + self.dicmap[pos_size]
    return  "%s_5p\t%s\n%s_3p\t%s" % (self.name, mir5p, self.name, mir3p)


mirdict = get_fasta (sys.argv[2])
dicobject = {}
for mir in mirdict:
  dicobject[mir] = Mirna(mir, mirdict[mir])
   
F = open (sys.argv[1], "r")
for line in F:
  fields = line.split()
  name = fields[1]
  offset= int(fields[2])
  sequence= fields[3]
  dicobject[name].addread(offset, len(sequence))
F.close()

F = open (sys.argv[4], "w")
for mir in sorted(dicobject):
  print >> F, dicobject[mir].hitmap()
  for i, counts in enumerate (dicobject[mir].density()):
    print >> F, "%s\t%s" % (i+1, counts) # attention 1-based offset for biologists
  print >> F
  print >> F, dicobject[mir].normalized_density()
F.close()
  
F = open (sys.argv[5], "w")
print >> F, "gene\t%s" % sys.argv[3]
for mir in sorted(dicobject):
  print >> F, dicobject[mir].splitcount_2(15)

