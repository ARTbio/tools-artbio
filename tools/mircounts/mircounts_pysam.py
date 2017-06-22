#!/usr/bin/env python

import pysam
import sys

bam_file = sys.argv[1]
gff_file = sys.argv[2]
print sys.argv 
samfile = pysam.AlignmentFile(bam_file, 'rb',check_sq=False)
headers = dict()
gff = open(gff_file,'r')
for line in gff.readlines():
    if line[0] != '#':
        gff_fields = line[:-1].split("\t")
        if gff_fields[2] == 'miRNA':
            headers[gff_fields[0]] = [samfile.count(reference=gff_fields[0], start=int(gff_fields[3]), end=int(gff_fields[4])), samfile.count_coverage(reference=gff_fields[0], start=int(gff_fields[3]), end=int(gff_fields[4]))]
gff.close()
for chrom in headers.keys():
    for pos in range(len(headers[chrom][1][1])):
        print "\t".join([chrom,str(headers[chrom][0]),str(pos),str(headers[chrom][1][0][pos]+headers[chrom][1][1][pos]+headers[chrom][1][2][pos]+headers[chrom][1][3][pos])])
