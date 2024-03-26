#!/usr/bin/env python
import argparse
import csv
import os
import shlex
import subprocess
import sys
from collections import defaultdict


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='''
             Prepartion of repetive element pseudogenomes bowtie\
             indexes and annotation files used by RepEnrich.py enrichment.''',
                                 prog='getargs_genome_maker.py')
parser.add_argument('--annotation_file', action='store',
                    metavar='annotation_file',
                    help='''Repeat masker annotation of the genome of\
                         interest. Download from RepeatMasker.org\
                         Example: mm9.fa.out''')
parser.add_argument('--genomefasta', action='store', metavar='genomefasta',
                    help='''Genome of interest in fasta format.\
                         Example: mm9.fa''')
parser.add_argument('--setup_folder', action='store', metavar='setup_folder',
                    help='''Folder that contains bowtie indexes of repeats and\
                         repeat element psuedogenomes.\
                         Example working/setup''')
parser.add_argument('--gaplength', action='store', dest='gaplength',
                    metavar='gaplength', default='200', type=int,
                    help='''Length of the N-spacer in the\
                         repeat pseudogenomes.  Default 200''')
parser.add_argument('--flankinglength', action='store', dest='flankinglength',
                    metavar='flankinglength', default='25', type=int,
                    help='''Length of the flanking regions used to build\
                         repeat pseudogenomes. Flanking length should be set\
                         according to the length of your reads.\
                         Default 25, for 50 nt reads''')
args = parser.parse_args()

# parameters from argsparse
gapl = args.gaplength
flankingl = args.flankinglength
annotation_file = args.annotation_file
genomefasta = args.genomefasta
setup_folder = args.setup_folder

# check that the programs we need are available
try:
    subprocess.call(shlex.split("bowtie --version"),
                    stdout=open(os.devnull, 'wb'),
                    stderr=open(os.devnull, 'wb'))
except OSError:
    print("Error: Bowtie not available in the path")
    raise


def starts_with_numerical(list):
    try:
        if len(list) == 0:
            return False
        return True
    except ValueError:
        return False


# Define a text importer
def import_text(filename, separator):
    csv.field_size_limit(sys.maxsize)
    file = csv.reader(open(filename), delimiter=separator,
                      skipinitialspace=True)
    return [line for line in file if starts_with_numerical(line)]


# Make a setup folder
if not os.path.exists(setup_folder):
    os.makedirs(setup_folder)

# load genome into dictionary and compute length
g = SeqIO.to_dict(SeqIO.parse(genomefasta, "fasta"))
genome = defaultdict(dict)

for chr in g.keys():
    genome[chr]['sequence'] = g[chr].seq
    genome[chr]['length'] = len(g[chr].seq)

# Build a bedfile of repeatcoordinates to use by RepEnrich region_sorter
repeat_elements = set()
rep_coords = defaultdict(list)  # Merged dictionary for coordinates

with open(os.path.join(setup_folder, 'repnames.bed'), 'w') as fout:
    f_in = import_text(annotation_file, ' ')
    for line in f_in:
        repname = line[9].translate(str.maketrans('()/', '___'))
        repeat_elements.add(repname)
        repchr, repstart, repend = line[4], line[5], line[6]
        fout.write(f"{repchr}\t{repstart}\t{repend}\t{repname}\n")
        rep_coords[repname].extend([repchr, repstart, repend])
# repeat_elements now contains the unique repeat names
# rep_coords is a dictionary where keys are repeat names and values are lists
# containing chromosome, start, and end coordinates for each repeat instance

# sort repeat_elements and print them in repgenomes_key.txt
with open(os.path.join(setup_folder, 'repgenomes_key.txt'), 'w') as fout:
    for i, repeat in enumerate(sorted(repeat_elements)):
        fout.write('\t'.join([repeat, str(i)]) + '\n')

# generate spacer for pseudogenomes
spacer = ''.join(['N' for i in range(gapl)])

# generate metagenomes and save them to FASTA files for bowtie build
for repname in rep_coords:
    metagenome = ''
    # iterating coordinate list by block of 3 (chr, start, end)
    block = 3
    for i in range(0, len(rep_coords[repname]) - block + 1, block):
        batch = rep_coords[repname][i:i+block]
        print(batch)
        chromosome = batch[0]
        start = max(int(batch[1]) - flankingl, 0)
        end = min(int(batch[2]) + flankingl,
                  int(genome[chromosome]['length'])-1) + 1
        metagenome = (
            f"{metagenome}{spacer}"
            f"{genome[chromosome]['sequence'][start:end]}"
            )

    # Create Fasta of repeat pseudogenome
    fastafilename = f"{os.path.join(setup_folder, repname)}.fa"
    record = SeqRecord(Seq(metagenome), id=repname, name='', description='')
    SeqIO.write(record, fastafilename, "fasta")

    # Generate repeat pseudogenome bowtie index
    bowtie_build_cmd = ["bowtie-build", "-f", fastafilename,
                        os.path.join(setup_folder, repname)]
    subprocess.run(bowtie_build_cmd, check=True)
