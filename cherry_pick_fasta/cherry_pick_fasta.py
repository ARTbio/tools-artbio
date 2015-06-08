#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Chery pick of fasta sequences satisfying a query string in their header/name
"""

import argparse

def Parser():
    the_parser = argparse.ArgumentParser(
        description="Cherry pick fasta sequences")
    the_parser.add_argument(
        '--input', action="store", type=str, help="input fasta file")
    the_parser.add_argument(
        '--query-string', dest="query_string", action="store", type=str,
                            help="header containing the string will be extracted as well as the corresponding sequence")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output fasta file")
    args = the_parser.parse_args()
    return args

def __main__():
    """ main function """
    args = Parser()
    search_term = args.query_string
    CrudeFasta = open (args.input, "r").read()
    Output = open (args.output, "w")
    FastaListe = CrudeFasta.split(">")
    for sequence in FastaListe:
        if search_term in sequence:
            print >> Output,  ">%s" % sequence.rstrip()
    Output.close()


if __name__ == "__main__":
    __main__()
