#!/usr/bin/python
# script for computing overlap signatures from a bowtie output
# Christophe Antoniewski <drosofff@gmail.com>
# Usage signature.py <1:input> <2:format of input> <3:minsize query> <4:maxsize query> <5:minsize target> <6:maxsize target>
#			  <7:minscope> <8:maxscope> <9:output> <10:bowtie index> <11:procedure option> <12: graph (global or lattice)>
# 			  <13: R code>
# version 2.0.0

import sys
import subprocess
import argparse
from smRtools import *
from collections import defaultdict  # test whether it is required


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="input alignment file")
    the_parser.add_argument('--inputFormat', action="store", type=str, choices=[
                            "tabular", "bam", "sam"], help="format of alignment file (tabular/bam/sam)")
    the_parser.add_argument(
        '--minquery', type=int, help="Minimum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxquery', type=int, help="Maximum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--mintarget', type=int, help="Minimum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxtarget', type=int, help="Maximum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--minscope', type=int, help="Minimum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--maxscope', type=int, help="Maximum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--outputOverlapDataframe', action="store", type=str, help="Overlap dataframe")
    the_parser.add_argument('--referenceGenome', action='store',
                            help="path to the bowtie-indexed or fasta reference")
    the_parser.add_argument('--extract_index', action='store_true',
                            help="specify if the reference is an indexed Bowtie reference")
    the_parser.add_argument('--graph', action='store', choices=[
                            "global", "lattice"], help="small RNA signature is computed either globally or by item (global-lattice)")
    the_parser.add_argument(
        '--rcode', type=str, help="R code to be passed to the python script")
    args = the_parser.parse_args()
    return args

args = Parser()

if args.extract_index:
    GenomeFormat = "bowtieIndex"
else:
    GenomeFormat = "fastaSource"

if args.inputFormat == "tabular":
    Genome = HandleSmRNAwindows(
        args.input, args.inputFormat, args.referenceGenome, GenomeFormat)
elif args.inputFormat == "sam":
    Genome = HandleSmRNAwindows(
        args.input, args.inputFormat, args.referenceGenome, GenomeFormat)
else:
    Genome = HandleSmRNAwindows(
        args.input, args.inputFormat, args.referenceGenome, GenomeFormat)

# replace objDic by Genome.instanceDict or... objDic = Genome.instanceDict
objDic = Genome.instanceDict

args.maxscope += 1

general_frequency_table = dict(
    [(i, 0) for i in range(args.minscope, args.maxscope)])
general_percent_table = dict(
    [(i, 0) for i in range(args.minscope, args.maxscope)])
OUT = open(args.outputOverlapDataframe, "w")

if args.graph == "global":
    # for normalized summing of local_percent_table(s)
    readcount_dic = {}
    Total_read_in_objDic = 0
    for item in objDic:
        readcount_dic[item] = objDic[item].readcount(
            args.minquery, args.maxquery)
        Total_read_in_objDic += readcount_dic[item]
    ######
    for x in (objDic):
        local_frequency_table = objDic[x].signature(
            args.minquery, args.maxquery, args.mintarget, args.maxtarget, range(args.minscope, args.maxscope))
        local_percent_table = objDic[x].hannon_signature(
            args.minquery, args.maxquery, args.mintarget, args.maxtarget, range(args.minscope, args.maxscope))
        try:
            for overlap in local_frequency_table.keys():
                general_frequency_table[overlap] = general_frequency_table.get(
                    overlap, 0) + local_frequency_table[overlap]
        except:
            pass
        try:
            for overlap in local_percent_table.keys():
                general_percent_table[overlap] = general_percent_table.get(
                    overlap, 0) + (1. / Total_read_in_objDic * readcount_dic[x] * local_percent_table[overlap])
        except:
            pass
    print >> OUT, "overlap\tnum of pairs\tprobability"
    for classe in sorted(general_frequency_table):
        print >> OUT, "%i\t%i\t%f" % (
            classe, general_frequency_table[classe], general_percent_table[classe])

else:
    print >> OUT, "overlap\tnum of pairs\tprobability\titem"
    for x in (objDic):
        local_frequency_table = objDic[x].signature(
            args.minquery, args.maxquery, args.mintarget, args.maxtarget, range(args.minscope, args.maxscope))
        local_percent_table = objDic[x].hannon_signature(
            args.minquery, args.maxquery, args.mintarget, args.maxtarget, range(args.minscope, args.maxscope))
        for classe in range(args.minscope, args.maxscope):
            print >> OUT, "%i\t%i\t%f\t%s" % (
                classe, local_frequency_table[classe], local_percent_table[classe], x)

OUT.close()

# Run the R script that is defined in the xml using the Rscript binary
# provided with R.
R_command = "Rscript " + args.rcode
process = subprocess.Popen(R_command.split())
