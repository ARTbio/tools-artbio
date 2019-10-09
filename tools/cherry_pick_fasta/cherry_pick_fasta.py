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
                            help="with, without, or withlist, withoutlist")
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


def parse_fasta_with(query, FastaListe):
    if not isinstance(query, list):
        query = [query]
    accumulator = []
    for sequence in FastaListe:
        for string in query:
            if string in sequence:
                accumulator.append(sequence)
                continue
    return accumulator


def complement_fasta(fullfasta, subfasta):
    return list(set(fullfasta) - set(subfasta))


def __main__():
    """ main function """
    args = Parser()
    search_term = args.query_string
    CrudeFasta = open(args.input, "r").read()
    Output = open(args.output, "w")
    FastaListe = CrudeFasta.split(">")[1:]
    if args.searchfor == 'with':
        contList = parse_fasta_with(search_term, FastaListe)
        contFasta = ">%s" % ">".join(contList)
        Output.write(contFasta)
    elif args.searchfor == 'without':
        notcontList = complement_fasta(FastaListe,
                                       parse_fasta_with(search_term,
                                                        FastaListe))
        notcontFasta = ">%s" % ">".join(notcontList)
        Output.write(notcontFasta)
    Output.close()


if __name__ == "__main__":
    __main__()
