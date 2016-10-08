#!/usr/bin/env python

import sys
import requests
import argparse
import csv


def get_runinfo(query):
    params = {'db': 'sra', 'term': query, 'usehistory': 'y', 'retmode': 'json'}
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    answer = requests.get(url, params=params).json()
    webenv = answer["esearchresult"]["webenv"]
    query_key = answer["esearchresult"]["querykey"]

    params = {'db': 'sra', 'save': 'efetch', 'rettype': 'runinfo', 'WebEnv': webenv, 'query_key': query_key}
    url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi'
    response = requests.get(url, params=params)
    return response.text.strip()


def join_response(response_list):
    """
    Joins responses from multiple esearch query responses
    """
    bodies = []
    for i, response in enumerate(response_list):
        lines = response.split('\n')
        if i == 0:
            header = [lines[0]]
        bodies.extend(lines[1:])
    return header + bodies


def write_tabular(joint_response, output):
    transformed_response = [line for line in csv.reader(joint_response)]
    [output.write("\t".join(line)+'\n') for line in transformed_response]


def process_query(queries, out_file):
    response_list = [get_runinfo(query) for query in queries]
    joint_response = join_response(response_list)
    if out_file == sys.stdout:
        output = sys.stdout
        write_tabular(joint_response, output)
    else:
        with open(out_file, 'w') as output:
            write_tabular(joint_response, output)


def parser():
    the_parser = argparse.ArgumentParser(
        description="Get runinfo for SRA or SRR accession")
    the_parser.add_argument(
        '--query', action="store", type=str, nargs='+',
        help="One or more SRA/SRR accessions, seperated by space", required=True)
    the_parser.add_argument(
        '--output', action="store", type=str, default=sys.stdout, help="output file")
    args = the_parser.parse_args()
    return args

args = parser()


def __main__():
    process_query(args.query, args.output)

if __name__ == "__main__":
    __main__()
