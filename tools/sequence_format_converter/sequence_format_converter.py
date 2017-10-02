#!/usr/bin/env python
#
import argparse
import logging
import sys
from collections import defaultdict


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str,
        help="input file, accepted format: fastq, fasta, fasta_weigthed, \
            tabular")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output converted file")
    the_parser.add_argument(
        '--format', action="store", type=str,
        help="select output format (fasta, fasta_weigthed, tabular")
    args = the_parser.parse_args()
    return args


class Sequencing:

    def __init__(self, input, output, format):
        self.input = input
        self.output = open(output, 'w')
        self.outputformat = format
        self.inputformat = self.detectformat(self.input)
        self.seqdic = defaultdict(int)
        self.read(self.input, self.inputformat)
        self.write(self.output, self.outputformat)

    def detectformat(self, input):
        input = open(input, 'r')
        block = []
        reference = ['A', 'T', 'G', 'C', 'N']
        format = ''
        try:
            for l in range(4):
                block.append(input.readline()[:-1])
        except:
            logging.info("File hasn't at leat four lines !")
            sys.exit("File hasn't at leat four lines !")
        input.close()
        line1, line2, line3, line4 = block[0], block[1], block[2], block[3]
        if line1[0] == '>' and line3[0] == '>':
            logging.info("'>' detected in lines 1 and 3")
            sequence = ''.join([line2, line4]).upper()
            nucleotides = set([base for base in sequence])
            for nucleotide in nucleotides:
                if nucleotide not in reference:
                    logging.info("But other nucleotides that A, T, G, C or N")
                    sys.exit('input appears to be Fasta but with \
                              unexpected nucleotides')
            format = 'fasta'
        elif line1[0] == '>' and line4[0] == '>':
            logging.info("'>' detected in lines 1 and 4")
            sequence = ''.join([line2, line3]).upper()
            nucleotides = set([base for base in sequence])
            for nucleotide in nucleotides:
                if nucleotide not in reference:
                    logging.info("But other nucleotides that A, T, G, C or N")
                    sys.exit('input appears to be Fasta but with \
                              unexpected nucleotides')
            format = 'fasta'
        elif line1[0] == '>':
            logging.info("'>' detected in lines 1")
            sequence = ''.join([line2, line3, line4]).upper()
            nucleotides = set([base for base in sequence])
            for nucleotide in nucleotides:
                if nucleotide not in reference:
                    logging.info("But other nucleotides that A, T, G, C or N")
                    sys.exit('input appears to be Fasta but with \
                              unexpected nucleotides')
            format = 'fasta'
        if format == 'fasta':
            try:
                for line in block:
                    if line[0] == '>':
                        int(line.split('_')[-1])
                return 'fastaw'
            except:
                return 'fasta'
        if line1[0] == '@' and line3[0] == '+':
            nucleotides = set([base for base in line2])
            for nucleotide in nucleotides:
                if nucleotide not in reference:
                    logging.info("Looks like fastq input but other nucleotides \
                                 that A, T, G, C or N")
                    sys.exit("input appears to be Fastq \
                             but with unexpected nucleotides")
            return 'fastq'
        for line in block:
            if len(line.split('\t')) != 2:
                logging.info("No valid format detected")
                sys.exit('No valid format detected')
            try:
                int(line.split('\t')[-1])
            except:
                logging.info("No valid format detected")
                sys.exit('No valid format detected')
            for nucleotide in line.split('\t')[0]:
                if nucleotide not in reference:
                    logging.info("No valid format detected")
                    sys.exit('No valid format detected')
        return 'tabular'

    def read(self, input, format):
        input = open(input, 'r')
        if format == 'fasta':
            try:
                self.readfasta(input)
            except:
                logging.info("an error occured while reading fasta")
        elif format == 'fastaw':
            try:
                self.readfastaw(input)
            except:
                logging.info("an error occured while reading fastaw")
        elif format == 'tabular':
            try:
                self.readtabular(input)
            except:
                logging.info("an error occured while reading tabular")
        elif format == 'fastq':
            try:
                self.readfastq(input)
            except:
                logging.info("an error occured while reading fastq")
        else:
            logging.info("no valid format detected")
            sys.exit('No valid format detected')

    def readfastaw(self, input):
        for line in input:
            if line[0] == ">":
                weigth = int(line[:-1].split("_")[-1])
            else:
                self.seqdic[line[:-1]] += weigth
        input.close()

    def readfasta(self, input):
        ''' this method is able to read multi-line fasta sequence'''
        for line in input:
            if line[0] == ">":
                try:
                    #  to dump the sequence of the previous item
                    #  try because of first missing stringlist variable
                    self.seqdic["".join(stringlist)] += 1
                except NameError:
                    pass
                stringlist = []
            else:
                try:
                    stringlist.append(line[:-1])
                except UnboundLocalError:
                    # if file went through filter and contains only empty lines
                    logging.info("first line is empty.")
        try:
            self.seqdic["".join(stringlist)] += 1  # for the last sequence
        except NameError:
            logging.info("input file has not fasta sequences.")
        input.close()

    def readtabular(self, input):
        for line in input:
            fields = line[:-1].split('\t')
            self.seqdic[fields[0]] += int(fields[1])
        input.close()

    def readfastq(self, input):
        linecount = 0
        for line in input:
            linecount += 1
            if linecount % 4 == 2:
                self.seqdic[line[:-1]] += 1
        input.close()

    def write(self, output, format='fasta'):
        if format == 'fasta':
            headercount = 0
            for seq in sorted(self.seqdic, key=self.seqdic.get, reverse=True):
                for i in range(self.seqdic[seq]):
                    headercount += 1
                    output.write('>%s\n%s\n' % (headercount, seq))
        elif format == 'fastaw':
            headercount = 0
            for seq in sorted(self.seqdic, key=self.seqdic.get, reverse=True):
                headercount += 1
                output.write('>%s_%s\n%s\n' % (headercount,
                                               self.seqdic[seq], seq))
        elif format == 'tabular':
            for seq in sorted(self.seqdic, key=self.seqdic.get, reverse=True):
                output.write('%s\t%s\n' % (seq, self.seqdic[seq]))
        output.close()


def main(input, output, format):
    Sequencing(input, output, format)


if __name__ == "__main__":
    args = Parser()
    log = logging.getLogger(__name__)
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    main(args.input, args.output, args.format)
