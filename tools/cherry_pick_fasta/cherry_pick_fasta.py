#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Chery pick of fasta sequences satisfying a query string in their header/name
"""

import argparse


def Parser():
    the_parser = argparse.ArgumentParser(
        description="Cherry pick fasta sequences")
    the_parser.add_argument('--input', action="store", type=str,
                            help="input fasta file")
    the_parser.add_argument('--searchfor', action="store", type=str,
                            help="with or without")
    the_parser.add_argument('--query-string', dest="query_string",
                            action="store", type=str,
                            help="headers containing the string will be \
                                  extracted or excluded as well as the \
                                  corresponding sequence")
    the_parser.add_argument('--query-file', dest="query_file",
                            action="store", type=str,
                            help="headers containing any of the strings provided in the \
                                  text file (1 string per line) will be \
                                  extracted or excluded as well as the \
                                  corresponding sequence")

    the_parser.add_argument(
        '--output', action="store", type=str, help="output fasta file")
    args = the_parser.parse_args()
    return args

def parse_fasta_with(query_string, FastaListe):
    accumulator = []
    for sequence in FastaListe:
        if query_string in sequence:
            accumulator.append(">%s\n" % sequence.rstrip())
    return "".join(accumulator)

def parse_fasta_without(query_string, FastaListe):
    accumulator = []
    for sequence in FastaListe:
        if query_string not in sequence:
            accumulator.append(">%s\n" % sequence.rstrip())
    return "".join(accumulator)


def __main__():
    """ main function """
    args = Parser()
    search_term = args.query_string
    CrudeFasta = open(args.input, "r").read()
    Output = open(args.output, "w")
    FastaListe = CrudeFasta.split(">")[1:]
    if args.searchfor == 'with':
        Output.write(parse_fasta_with(search_term, FastaListe))
    else:
        Output.write(parse_fasta_without(search_term, FastaListe))
    Output.close()


if __name__ == "__main__":
    __main__()
