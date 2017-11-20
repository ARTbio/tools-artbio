#!/usr/bin/env python

from subprocess import call
import re
import argparse
import sys
import os


def command_parse():
    LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    parser = argparse.ArgumentParser(description='fetch_fasta_from_NCBI runner')
    parser.add_argument('-i', dest='query_string', help='NCBI Query String',
                        required=True)
    parser.add_argument('-t', dest='tool_path', help='fetch_fasta_from_NCBI.py\
    path', required=True)
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

def get_number_of_uids(args):
    '''
    This function launches a dry_run and looks for the number of UIDs in the log
    '''
    tmp_logfile = 'tmp_logfile.log'
    tmp_stderrfile = 'tmp_stderrfile.err'
    pattern = re.compile(".*Found\s+(\d+)\s+UID.*")
    return_code = 0
    number_of_uids = 0
    try:
        # Launch dry_run
        with open(tmp_stderrfile, 'w') as ferr:
            dry_run_cmd = ['python', args.tool_path, '-i',
                           args.query_string, '-c', '-l', tmp_logfile,
                           '-d', args.dbname]
            return_code = call(dry_run_cmd, stderr=ferr)
        assert return_code == 0 # Check if return code is OK
        with open(tmp_logfile, 'r') as f:# Get UID count
            for line in f:
                if pattern.match(line):
                    number_of_uids = pattern.match(line).group(1)
        return number_of_uids
    except AssertionError as e:# Returncode is not 0
        try:# Check if the returncode isn't 0 because no UIDs were found
            with open(tmp_logfile, 'r') as f:
                for line in f:
                    if pattern.match(line):
                        number_of_uids = pattern.match(line).group(1)
            os.remove(tmp_logfile)
            return number_of_uids
        except IOError as e:# No log file was found so there is a real problem
            to_stderr = ""
            with open(tmp_stderrfile, 'w') as f:
                to_stderr = f.read()
        os.remove(tmp_stderrfile)
        raise ValueError('Non zero returncode: %s\n%s', return_code, to_stderr)

def __main__():
    args = command_parse()
    print(str(get_number_of_uids(args)))

if __name__ == '__main__':
    __main__()
