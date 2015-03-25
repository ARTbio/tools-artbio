#!/usr/bin/python
# yac = yet another clipper
# v 1.1.0 - 23-08-2014 - argparse implementation
# Usage yac.py  $input $output $adapter_to_clip $min $max $Nmode
# Christophe Antoniewski <drosofff@gmail.com>

import sys, string, argparse

def Parser():
  the_parser = argparse.ArgumentParser()
  the_parser.add_argument('--input', action="store", nargs='+', help="input fastq files")
  the_parser.add_argument('--output', action="store", type=str, help="output, clipped fasta file")
  the_parser.add_argument('--adapter_to_clip', action="store", type=str, help="adapter sequence to clip")
  the_parser.add_argument('--min', action="store", type=int, help="minimal size of clipped sequence to keep")
  the_parser.add_argument('--max', action="store", type=int, help="maximal size of clipped sequence to keep")
  the_parser.add_argument('--Nmode', action="store", type=str, choices=["accept", "reject"], help="accept or reject sequences with N for clipping")
  args = the_parser.parse_args()
  args.adapter_to_clip = args.adapter_to_clip.upper()
  return args


class Clip:
  def __init__(self, inputfile, outputfile, adapter, minsize, maxsize):
    self.inputfile = inputfile
    self.outputfile = outputfile
    self.adapter = adapter
    self.minsize = int(minsize)
    self.maxsize = int(maxsize)
    def motives (sequence):
      '''return a list of motives for perfect (6nt) or imperfect (7nt with one mismatch) search on import string module'''
      sequencevariants = [sequence[0:6]] # initializes the list with the 6mer perfect match
      dicsubst= {"A":"TGCN", "T":"AGCN", "G":"TACN", "C":"GATN"}
      for pos in enumerate(sequence[:6]):
        for subst in dicsubst[pos[1]]:
          sequencevariants.append(sequence[:pos[0]]+ subst + sequence[pos[0]+1:7])
      return sequencevariants
    self.adaptmotifs= motives(self.adapter)

  def scanadapt(self, adaptmotives=[], sequence=""):
    '''scans sequence for adapter motives'''
    if sequence.rfind(adaptmotives[0]) != -1:
      return sequence[:sequence.rfind(adaptmotives[0])]
    for motif in adaptmotives[1:]:
      if sequence.rfind(motif) != -1:
        return sequence[:sequence.rfind(motif)]
    return sequence

  def clip_with_N (self):
    '''clips adapter sequences from inputfile. 
    Reads containing N are retained.'''
    global id
    iterator = 0
    F = open (self.inputfile, "r")
    O = open (self.outputfile, "a")
    for line in F:
      iterator += 1
      if iterator % 4 == 2:
        trim = self.scanadapt (self.adaptmotifs, line.rstrip() )
        if self.minsize <= len(trim) <= self.maxsize:
          id += 1
          print >> O, ">%i\n%s" % (id, trim)
    F.close()
    O.close()
  def clip_without_N (self):
    '''clips adapter sequences from inputfile. 
    Reads containing N are rejected.'''
    global id
    iterator = 0
    F = open (self.inputfile, "r")
    O = open (self.outputfile, "a")
    for line in F:
      iterator += 1
      if iterator % 4 == 2:
        trim = self.scanadapt (self.adaptmotifs, line.rstrip() )
        if "N" in trim: continue
        if self.minsize <= len(trim) <= self.maxsize:
          id += 1
          print >> O, ">%i\n%s" % (id, trim)
    F.close()
    O.close()

def __main__ (inputfile, outputfile, adapter, minsize, maxsize, Nmode):
  instanceClip = Clip (inputfile, outputfile, adapter, minsize, maxsize)
  if Nmode == "accept":
    instanceClip.clip_with_N()
  else:
    instanceClip.clip_without_N()

if __name__ == "__main__" :
  args = Parser()
  id = 0
  for inputfile in args.input:
    __main__(inputfile, args.output, args.adapter_to_clip, args.min, args.max, args.Nmode)
