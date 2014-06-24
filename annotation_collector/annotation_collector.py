#!/usr/bin/env python
#By, drosofff@gmail.com
# command: annotation_collector.py $output input1, label1, input2, label2, etc...

import sys, os

def countlineinfile(file):
  F = open (file, "r")
  count = 0
  for line in F:
    count += 1
  F.close()
  return count/2
results = []

for file, label in zip (sys.argv[2:-1:2], sys.argv[3:-1:2]):
  results.append ( (countlineinfile(file), label) )
  
Fout = open (sys.argv[1], "w")

print >> Fout, "# %s" % (sys.argv[-1])
for filecount, label in results:
  print >> Fout, "%s\t%s\t%.2f %%" % (label, filecount, filecount/float(results[0][0])*100 )
print >> Fout  
Fout.close()
