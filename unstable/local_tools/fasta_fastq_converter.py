#!/usr/bin/python
#
import sys

input = open(sys.argv[1], "r")
output = open(sys.argv[2], "w")

for line in input:
  if line[0]==">":
    print >> output, "@HTW-"+line[1:-1]
    continue
  else:
    print >> output, line[:-1]
    print >> output, "+"
    print >> output, "H"*len(line[:-1])
    
input.close()
output.close()
