#!/usr/bin/env python3

import optparse


def parse_options():
    """
    Parse the options guiven to the script
    """
    parser = optparse.OptionParser(description='Get unmatched blast queries')
    parser.add_option('-f','--fasta', dest='fasta_file', help='Query fasta file\
used during blast')
    parser.add_option('-b','--blast', dest='blast_file', help='Blast tabular\
output (queries in 1rst column)')
    parser.add_option('-o','--output', dest='output_file', help='Output file name')
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')
    return options

def get_matched(blast_file):
    """
    Get a dictionary of all the queries that got a match
    """
    matched = dict()
    blast_file_handle = open(blast_file, 'r')
    for line in blast_file_handle.readlines():
        fields = line.split("\t")
        query_id = fields[0]
        matched[query_id] = 1
    blast_file_handle.close()
    return matched

def get_unmatched(output_file, fasta_file, matched):
    """
    Compares matched queries to query fasta file and print unmatched to ouput
    """
    output_file_handle = open(output_file, 'w')
    fasta_file_handle = open(fasta_file, 'r')
    unmatched = False
    for line in fasta_file_handle.readlines():
        if line.startswith('>'):
            subline = line[1:100].rstrip() #qid are 100chars long in blast
            if subline not in matched:
                output_file_handle.write(line)
                unmatched = True
            else:
                unmatched = False
        elif unmatched:
            output_file_handle.write(line)
    fasta_file_handle.close()
    output_file_handle.close()

def __main__():
    opts = parse_options()
    matched = get_matched(opts.blast_file)
    get_unmatched(opts.output_file, opts.fasta_file, matched)

if __main__():
    __main__()
