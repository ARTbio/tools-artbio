#!/usr/bin/python
# script preparing tab file to plot in R for small RNA profiling
# version 1 29-1-2012
# Usage plotter.py <bowtie input> <min size> <max size> <normalization factor> <tabular output>


import sys

def acquisition (file2parse, sizerange):
  F = open (file2parse)
  plus_table = {}
  minus_table = {}
  for line in F:
    field = line.split()
    coordinate = int( field[3] )
    strand =  field[1]
    sequence = field[4]
    size = len (sequence )
    if strand == "+" and size in sizerange:
      plus_table[coordinate] = plus_table.get(coordinate, 0) + 1
    if strand == "-" and size in sizerange:
      coordinate = coordinate + size -1 # 23-11-2012 : this line was missing ! it is a BUG that probably altered the Nature maps :-((
      minus_table[coordinate] = minus_table.get(coordinate, 0) + 1
  return plus_table, minus_table
    
  
def output_table (plus_table, minus_table, Nfactor, output):
  Nfactor = float(Nfactor)
  plus_coordinates = set( plus_table.keys() ) 
  minus_coordinates = set( minus_table.keys() ) 
  coords = sorted (plus_coordinates.union (minus_coordinates) )
  ## added 23-2-2013 to have, instead, exaustive coordinates
  ## coords = range (min(coords), max(coords) + 1)
  ##
  OUT = open (output, "w")
  print >> OUT, "coord\tplus\tminus"
  for coordinate in coords :
    print >> OUT, "%s\t%s\t%s" % ( coordinate, plus_table.get(coordinate, 0)*Nfactor, - minus_table.get(coordinate, 0)*Nfactor )
  

def sizing (minsize, maxsize) :
  size_range = range ( int (minsize), int (maxsize) + 1 )
  return size_range
  
plus_table, minus_table = acquisition (sys.argv[1], sizing ( sys.argv[2], sys.argv[3] ) )
output_table ( plus_table, minus_table, sys.argv[4], sys.argv[5] )
