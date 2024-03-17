#!/usr/bin/env python
import argparse
import csv
import os
import shlex
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='''
             Part I: Prepartion of repetive element pseudogenomes bowtie\
             indexes and annotation files used by RepEnrich.py\
             enrichment. Bowtie must be accessible and is loaded by the\
             script.  Repeat element pseudogenomes are indexed in order\
             to analyze reads that map to multiple locations of the genome.\
             The repeat element bamfiles are prepared in order to use a\
             region sorter to analyze reads that map to a single location\
             of the genome.''',
                                 prog='getargs_genome_maker.py')
parser.add_argument('--annotation_file', action='store',
                    metavar='annotation_file',
                    help='''This annotation file contains\
                         the repeat masker annotation for the genome of\
                         interest and may be downloaded at RepeatMasker.org\
                         Example: /data/annotation/mm9/mm9.fa.out''')
parser.add_argument('--genomefasta', action='store', metavar='genomefasta',
                    help='''Genome of interest in fasta format.\
                         Example: /data/annotation/mm9/mm9.fa''')
parser.add_argument('--setup_folder', action='store', metavar='setup_folder',
                    help='''Folder that contains bamfiles for repeats and\
                         repeat element psuedogenomes.\
                         Example /data/annotation/mm9/setup''')
parser.add_argument('--gaplength', action='store', dest='gaplength',
                    metavar='gaplength', default='200', type=int,
                    help='''Length of the spacer used to build\
                         repeat pseudogenomes.  Default 200''')
parser.add_argument('--flankinglength', action='store', dest='flankinglength',
                    metavar='flankinglength', default='25', type=int,
                    help='''Length of the flanking regions\
                         used to build repeat pseudogenomes.\
                         The flanking length should be set according to the\
                         length of your reads.  Default 25''')
args = parser.parse_args()

# parameters and paths specified in args_parse
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
    print("Error: Bowtie not loaded")
    raise

# Define a text importer
csv.field_size_limit(sys.maxsize)


def import_text(filename, separator):
    for line in csv.reader(open(os.path.realpath(filename)),
                           delimiter=separator, skipinitialspace=True):
        if line:
            yield line


# Make a setup folder
if not os.path.exists(setup_folder):
    os.makedirs(setup_folder)
# load genome into dictionary
print("loading genome...")
g = SeqIO.to_dict(SeqIO.parse(genomefasta, "fasta"))

print("Precomputing length of all chromosomes...")
idxgenome, lgenome, genome = {}, {}, {}
allchrs = g.keys()
for k, chr in enumerate(allchrs):
    genome[chr] = g[chr].seq
    lgenome[chr] = len(genome[chr])
    idxgenome[chr] = k

# Build a bedfile of repeatcoordinates to use by RepEnrich region_sorter
repeat_elements = []
# these dictionaries will contain lists
rep_chr, rep_start, rep_end = {}, {}, {}
fin = import_text(annotation_file, ' ')
with open(os.path.join(setup_folder, 'repnames.bed'), 'w') as fout:
    for i in range(3):
        next(fin)
    for line in fin:
        repname = line[9].translate(str.maketrans('()/', '___'))
        if repname not in repeat_elements:
            repeat_elements.append(repname)
        repchr = line[4]
        repstart = line[5]
        repend = line[6]
        fout.write('\t'.join([repchr, repstart, repend, repname]) + '\n')
        if repname in rep_chr:
            rep_chr[repname].append(repchr)
            rep_start[repname].append(repstart)
            rep_end[repname].append(repend)
        else:
            rep_chr[repname] = [repchr]
            rep_start[repname] = [repstart]
            rep_end[repname] = [repend]

# sort repeat_elements and print them in repgenomes_key.txt
with open(os.path.join(setup_folder, 'repgenomes_key.txt'), 'w') as fout:
    for i, repeat in enumerate(sorted(repeat_elements)):
        fout.write('\t'.join([repeat, str(i)]) + '\n')

# generate spacer for pseudogenomes
spacer = ''.join(['N' for i in range(gapl)])

# generate metagenomes and save them to FASTA files for bowtie build
for repname in rep_chr:
    metagenome = ''
    for i, repeat in enumerate(rep_chr[repname]):
        try:
            chromosome = rep_chr[repname][i]
            start = max(int(rep_start[repname][i]) - flankingl, 0)
            end = min(int(rep_end[repname][i]) + flankingl,
                      int(lgenome[chr])-1) + 1
            metagenome = f"{metagenome}{spacer}{genome[chromosome][start:end]}"
        except KeyError:
            print("Unrecognised Chromosome: " + rep_chr[repname][i])

    # Create the SeqRecord object
    record = SeqRecord(Seq(metagenome), id=repname, name='', description='')

    # Build the FASTA filename using f-string formatting
    fastafilename = f"{os.path.join(setup_folder, repname)}.fa"

    # Write the record to FASTA file directly using SeqIO.write()
    SeqIO.write(record, fastafilename, "fasta")

    # Construct the bowtie-build command as a list
    bowtie_build_cmd = ["bowtie-build", "-f",
                        fastafilename,
                        os.path.join(setup_folder, repname)]

    # Execute the command using subprocess.run()
    subprocess.run(bowtie_build_cmd, check=True)
