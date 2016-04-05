#!/usr/bin/env python
#
import argparse
import logging
import sys
from collections import defaultdict


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="input file")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output converted file")
    the_parser.add_argument(
        '--type', action="store", type=str, help="type of convertion")
    args = the_parser.parse_args()
    return args


def readfasta_writetabular(fasta, tabular, mode="oneline"):
    for line in fasta:
        if line[0] == ">":
            try:
                seqdic["".join(stringlist)] += 1 # to dump the sequence of the previous item - try because of first missing stringlist variable
            except NameError: pass
            stringlist=[]
        else:
            try:
                stringlist.append(line[:-1])
            except UnboundLocalError:  # if file went through filter and contains only empty lines
                logging.info("first line is empty.")
    try:
        seqdic["".join(stringlist)] +=  1 # for the last sequence
    except NameError:
        logging.info("input file has not fasta sequences.")
    for seq in sorted(seqdic, key=seqdic.get, reverse=True):
        tabular.write( "%s\t%s\n" % (seq, seqdic[seq]))


def readtabular_writefasta(tabular, fasta):
    counter = 0
    for line in tabular:
        fields = line.split()
        for i in range(int(fields[1])):
            counter += 1
            fasta.write( ">%s\n%s\n" % (counter, fields[0]) )


def readtabular_writefastaweighted (tabular, fasta):
    counter = 0
    for line in tabular:
        counter += 1
        fields = line[:-1].split()
        fasta.write( ">%s_%s\n%s\n" % (counter, fields[1],  fields[0]) )


def readfastaweighted_writefastaweighted(fastaweigthed_input, fastaweigthed_reparsed):
    number_reads = 0
    for line in fastaweigthed_input:
        if line[0] == ">":
            weigth = int(line[1:-1].split("_")[-1])
            number_reads += weigth
        else:
            seqdic[line[:-1]] += weigth
    n=0
    for seq in sorted(seqdic, key=seqdic.get, reverse=True):
        n += 1
        fastaweigthed_reparsed.write( ">%s_%s\n%s\n" % (n, seqdic[seq], seq) )
    log.info( "%s reads collapsed" % number_reads)


def readfastaweighted_writefasta(fastaweigthed, fasta):
    counter = 0
    for line in fastaweigthed:
        if line[0] == ">":
            weigth = int(line[1:-1].split("_")[-1])
        else:
            seq = line[:-1]
            for i in range (weigth):
                counter += 1
                fasta.write( ">%s\n%s\n" % (counter, seq) )


def main(input, output, type):
    with open(input, "r") as input:
        with open(output, "w") as output:
            if type == "fasta2tabular":
                readfasta_writetabular(input, output)
            elif type == "tabular2fasta":
                readtabular_writefasta(input, output)
            elif type == "tabular2fastaweight":
                readtabular_writefastaweighted (input, output)
            elif type == "fastaweight2fastaweight":
                readfastaweighted_writefastaweighted(input, output)
            elif type == "fastaweight2fasta":
                readfastaweighted_writefasta(input, output)


if __name__ == "__main__":
    seqdic = defaultdict(int)
    args = Parser()
    log = logging.getLogger(__name__)
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    main (args.input, args.output, args.type)
