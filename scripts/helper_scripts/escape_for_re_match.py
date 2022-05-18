#!/usr/bin/env python

## Use this script to generate escaped output if using re.match
## for verifying test-data. See deseq2 wrapper for output.

import sys

input = open(sys.argv[1])
output = open(sys.argv[2], 'w')

for i,line in enumerate(input):
    if i==0:
        output.write(line)
    else:
        newline = []
        fields=line.strip().split('\t')
        for i, field in enumerate(fields):
            if not i==0:
                if len(field)>3:
                    if 'e' in field:
                        field=field.split('e')[0][:-2]
                    field=field[:-3]+'.*?'
                    if '.' in field:
                        fields=field.split('\.')
                        field=fields[0]+'\.'+fields[1][:-5]+'.*?'
                else:
                    field='.*?'
            newline.append(field)
        output.write("\t".join(newline)+"\n")
input.close()
output.close()
