#!/usr/bin/python
# yac = yet another clipper
# v 1.2.1 - 23-08-2014 - Support FastQ output
# v 1.1.0 - 23-08-2014 - argparse implementation
# Usage yac.py  $input $output $adapter_to_clip $min $max $Nmode
# Christophe Antoniewski <drosofff@gmail.com>

import sys
import string
import argparse
from itertools import islice


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", nargs='+', help="input fastq files")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output, clipped fasta file")
    the_parser.add_argument(
        '--output_format', action="store", type=str, help="output format, fasta or fastq")
    the_parser.add_argument(
        '--adapter_to_clip', action="store", type=str, help="adapter sequence to clip")
    the_parser.add_argument(
        '--min', action="store", type=int, help="minimal size of clipped sequence to keep")
    the_parser.add_argument(
        '--max', action="store", type=int, help="maximal size of clipped sequence to keep")
    the_parser.add_argument('--Nmode', action="store", type=str, choices=[
                            "accept", "reject"], help="accept or reject sequences with N for clipping")
    args = the_parser.parse_args()
    args.adapter_to_clip = args.adapter_to_clip.upper()
    return args


class Clip:

    def __init__(self, inputfile, outputfile, output_format, adapter, minsize, maxsize, Nmode):
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.output_format = output_format
        self.adapter = adapter
        self.minsize = int(minsize)
        self.maxsize = int(maxsize)
        self.Nmode = Nmode

        def motives(sequence):
            '''return a list of motives for perfect (6nt) or imperfect (7nt with one mismatch) search on import string module'''
            sequencevariants = [
                sequence[0:6]]  # initializes the list with the 6mer perfect match
            dicsubst = {"A": "TGCN", "T": "AGCN", "G": "TACN", "C": "GATN"}
            for pos in enumerate(sequence[:6]):
                for subst in dicsubst[pos[1]]:
                    sequencevariants.append(
                        sequence[:pos[0]] + subst + sequence[pos[0] + 1:7])
            return sequencevariants
        self.adaptmotifs = motives(self.adapter)

    def scanadapt(self, adaptmotives=[], sequence="", qscore=""):
        '''scans sequence for adapter motives'''
        match_position = sequence.rfind(adaptmotives[0])
        if match_position != -1:
            return sequence[:match_position], qscore[:match_position]
        for motif in adaptmotives[1:]:
            match_position = sequence.rfind(motif)
            if match_position != -1:
                return sequence[:match_position], qscore[:match_position]
        return sequence, qscore

    def write_output(self, id, read, qscore, output):
        if self.output_format == "fasta":
            block = ">{0}\n{1}\n".format(id, read)
        else:
            block = "@HWI-{0}\n{1}\n+\n{2}\n".format(id, read, qscore)
        output.write(block)

    def handle_io(self):
        '''Open input file, pass read sequence and read qscore to clipping function.
        Pass clipped read and qscore to output function.'''
        id = 0
        output = open(self.outputfile, "a")
        with open(self.inputfile, "r") as input:
            block_gen = islice(input, 1, None, 2)
            for i, line in enumerate(block_gen):
                if i % 2:
                    qscore = line.rstrip()
                else:
                    read = line.rstrip()
                    continue
                trimmed_read, trimmed_qscore = self.scanadapt(
                    self.adaptmotifs, read, qscore)
                if self.minsize <= len(trimmed_read) <= self.maxsize:
                    if (self.Nmode == "reject") and ("N" in trimmed_read):
                        continue
                    id += 1
                    self.write_output(id, trimmed_read, trimmed_qscore, output)
            output.close()


def main(*argv):
    instanceClip = Clip(*argv)
    instanceClip.handle_io()

if __name__ == "__main__":
    args = Parser()
    id = 0
    for inputfile in args.input:
        main(inputfile, args.output, args.output_format,
             args.adapter_to_clip, args.min, args.max, args.Nmode)
