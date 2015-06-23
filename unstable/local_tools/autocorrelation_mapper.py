#!/usr/bin/python
# script to calculate autocorrelation profile of a window within in larger profile
# version 1 - 21-03-2013
# Usage autocorrelation_mapper.py <readmap> <window_up_coordinates> <window_down_coordinates> <output>

import sys
from scipy import stats

window_up_coor= int(sys.argv[2])
window_down_coor= int(sys.argv[3])

F = open (sys.argv[1], "r") # F is tabulated file coord \t forward \t reverse
lines = F.readlines()
start_coordinate = int (lines[2].split()[0]) # extract first coordinate from the second line, first field and ensure integer type.
end_coordinate = int (lines[-1].split()[0]) # extract last coordinate from the last line, first field and ensure interger type
reference_forward = []
reference_reverse = []
coordinates = []

pattern_forward =[]
pattern_reverse = []

for line in lines[1:]: # to skip the header line
  fields  = line.split()
  coordinates.append(int(fields[0]))
  reference_forward.append(int(float(fields[1])))
  reference_reverse.append(-int(float(fields[2])))
  if window_up_coor <= int(fields[0]) <= window_down_coor:
    pattern_forward.append(int(float(fields[1])))
    pattern_reverse.append(-int(float(fields[2])))
F.close()

# optional to find pattern on sense and antisense strand
pattern_reverse = pattern_forward[::-1]

OUT = open (sys.argv[4], "w")
 
windowsize = len(pattern_forward)

for coordinate in range(len(reference_forward) - windowsize):
  local_forward=reference_forward[coordinate:coordinate + windowsize]
  local_reverse=reference_reverse[coordinate:coordinate + windowsize]
  if sum(local_forward) + sum(local_reverse) < 50:
    continue
  if sum(local_forward) == 0:
    reference_to_local_cor_forward = [0,0]
  else:
    reference_to_local_cor_forward = stats.spearmanr(local_forward, pattern_forward)
  if sum(local_reverse) == 0:
    reference_to_local_cor_reverse = [0,0]
  else:
    reference_to_local_cor_reverse = stats.spearmanr(local_reverse, pattern_reverse)
  if (reference_to_local_cor_forward[0] > 0.2 or  reference_to_local_cor_reverse[0]>0.2): # here to find traffic jam like patterns on both sense and antisense of 2L
    print >> OUT, "%s\t%s\t%s\t%s" % (coordinates[coordinate], reference_to_local_cor_forward[0], reference_to_local_cor_reverse[0], sum(local_forward)+sum(local_reverse))
    OUT.flush()
