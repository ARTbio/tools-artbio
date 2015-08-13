#!/usr/bin/env python
# python parser module for pre-mir and mature miRNAs, guided by mirbase.org GFF3
# version 0.0.9 (1-6-2014)
# Usage MirParser.py  <1:index source> <2:extraction directive> <3:output pre-mir> <4: output mature miRs> <5:mirbase GFF3>
#                     <6:pathToLatticeDataframe or "dummy_dataframe_path"> <7:Rcode or "dummy_plotCode"> <8:latticePDF or "dummy_latticePDF">
#                     <9:10:11 filePath:FileExt:FileLabel> <.. ad  lib>

import sys
import subprocess

from smRtools import *

IndexSource = sys.argv[1]
ExtractionDirective = sys.argv[2]
if ExtractionDirective == "--do_not_extract_index":
    genomeRefFormat = "fastaSource"
elif ExtractionDirective == "--extract_index":
    genomeRefFormat = "bowtieIndex"
OutputPre_mirs = sys.argv[3]
OutputMature_Mirs = sys.argv[4]
GFF3_file = sys.argv[5]
lattice = sys.argv[6]
Rcode = sys.argv[7]
latticePDF = sys.argv[8]
Triplets = [sys.argv[9:][i:i + 3] for i in xrange(0, len(sys.argv[9:]), 3)]
MasterListOfGenomes = {}

for [filePath, FileExt, FileLabel] in Triplets:
    print FileLabel
    MasterListOfGenomes[FileLabel] = HandleSmRNAwindows(alignmentFile=filePath,
                                                        alignmentFileFormat=FileExt,
                                                        genomeRefFile=IndexSource,
                                                        genomeRefFormat=genomeRefFormat,
                                                        biosample=FileLabel)

header = ["gene"]
for [filePath, FileExt, FileLabel] in Triplets:
    header.append(FileLabel)

hit_table = ["\t".join(header)]  # table header: gene, sample1, sample2, sample3, etc. separated by tabulation

# read GFF3 to subinstantiate
gff3 = open(GFF3_file, "r")
lattice_dataframe = []
for line in gff3:
    if line[0] == "#":
        continue
    gff_fields = line[:-1].split("\t")
    chrom = gff_fields[0]
    gff_name = gff_fields[-1].split("Name=")[-1].split(";")[0]  # to isolate the GFF Name
    item_upstream_coordinate = int(gff_fields[3])
    item_downstream_coordinate = int(gff_fields[4])
    if gff_fields[6] == "+":
        item_polarity = "forward"
    else:
        item_polarity = "reverse"
    item_line = [gff_name]
    for sample in header[1:]:
        count = MasterListOfGenomes[sample].instanceDict[chrom].readcount(upstream_coord=item_upstream_coordinate,
                                                                          downstream_coord=item_downstream_coordinate,
                                                                          polarity=item_polarity)
        item_line.append(str(count))
        # subtreatement for lattice
        if lattice != "dummy_dataframe_path":
            if ("5p" not in gff_name) and ("3p" not in gff_name):
                lattice_dataframe.append(MasterListOfGenomes[sample].instanceDict[chrom].readcoverage(
                    upstream_coord=item_upstream_coordinate,
                    downstream_coord=item_downstream_coordinate,
                    windowName=gff_name + "_" + sample))
        # end of subtreatement for lattice
    hit_table.append("\t".join(item_line))
gff3.close()

Fpremirs = open(OutputPre_mirs, "w")
print >> Fpremirs, hit_table[0]
finalPreList = [i for i in sorted(hit_table[1:]) if ("5p" not in i) and ("3p" not in i)]
print >> Fpremirs, "\n".join(finalPreList)
Fpremirs.close()

Fmaturemires = open(OutputMature_Mirs, "w")
print >> Fmaturemires, hit_table[0]
finalMatureList = [i for i in sorted(hit_table[1:]) if ("5p" in i) or ("3p" in i)]
print >> Fmaturemires, "\n".join(finalMatureList)
Fmaturemires.close()

if lattice != "dummy_dataframe_path":
    Flattice = open(lattice, "w")
    print >> Flattice, "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("sample",
                                                       "mir",
                                                       "offset",
                                                       "offsetNorm",
                                                       "counts",
                                                       "countsNorm",
                                                       "polarity")
    print >> Flattice, "\n".join(lattice_dataframe)
    Flattice.close()
    R_command = "Rscript " + Rcode
    process = subprocess.Popen(R_command.split())
    process.wait()
