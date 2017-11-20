#!/usr/bin/env python

import sys
import argparse


def command_parse():
    LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    parser = argparse.ArgumentParser(description='fetch_fasta_from_NCBI runner')
    parser.add_argument('-i', dest='query_string', help='NCBI Query String',
                        required=True)
    parser.add_argument('-t', dest='tool_path', help='fetch_fasta_from_NCBI.py\
    path')
    parser.add_argument('-o', dest='outname', help='output file name')
    parser.add_argument('-d', dest='dbname', help='database type')
    parser.add_argument('--count', '-c', dest='count_ids',
                        action='store_true', default=False,
                        help='dry run ouputing only the number of sequences\
                        found')
    parser.add_argument('-l', '--logfile', help='log file (default=stderr)')
    parser.add_argument('--loglevel', choices=LOG_LEVELS, default='INFO',
                        help='logging level (default: INFO)')
    args = parser.parse_args()
    return args

def __main__():
    args = command_parse()
if tester:
    __main__()
