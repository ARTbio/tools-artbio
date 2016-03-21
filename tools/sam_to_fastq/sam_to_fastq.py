#!/usr/bin/python
#
import sys
import argparse

def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="input SAM file")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output FASTQ file")
    args = the_parser.parse_args()
    return args
    
        
def print_fastq_sequence(samline, file):
  samfields = samline[:-1].split("\t")
  file.write ( '@%s\n%s\n+\n%s\n' % (samfields[0], samfields[9], samfields[10]) )

def main(input, output):
    infile = open (input, "r")
    outfile = open (output, "w")
    with open (input, "r") as infile:
        with open (output, "w") as outfile:
            for line in infile:
                if line[0] == "@":
                    continue
                if line.split("\t")[1] != "4":
                    print_fastq_sequence (line, outfile)

if __name__ == "__main__":
    args = Parser()
    main (args.input, args.output)
