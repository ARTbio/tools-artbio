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
             Part I: Prepartion of repetive element psuedogenomes and repetive\
             element bamfiles.  This script prepares the annotation used by\
             downstream applications to analyze for repetitive element\
             enrichment. For this script to run properly bowtie must be\
             loaded.  The repeat element psuedogenomes are prepared in order\
             to analyze reads that map to multiple locations of the genome.\
             The repeat element bamfiles are prepared in order to use a\
             region sorter to analyze reads that map to a single location\
             of the genome. You will 1) annotation_file:\
             The repetitive element annotation file downloaded from\
             RepeatMasker.org database for your organism of interest.\
             2) genomefasta: Your genome of interest in fasta format,\
             3)setup_folder: a folder to contain repeat element setup files\
             command-line usage
             EXAMPLE: python master_setup.py\
             /users/nneretti/data/annotation/mm9/mm9_repeatmasker.txt\
             /users/nneretti/data/annotation/mm9/mm9.fa\
             /users/nneretti/data/annotation/mm9/setup_folder''',
                                 prog='getargs_genome_maker.py')
parser.add_argument('--annotation_file', action='store',
                    metavar='annotation_file',
                    help='''This annotation file contains\
                         the repeat masker annotation for the genome of\
                         interest and may be downloaded at RepeatMasker.org\
                         Example /data/annotation/mm9/mm9.fa.out''')
parser.add_argument('--genomefasta', action='store', metavar='genomefasta',
                    help='''Genome of interest in fasta format.\
                         Example /data/annotation/mm9/mm9.fa''')
parser.add_argument('--setup_folder', action='store', metavar='setup_folder',
                    help='''Folder that contains bamfiles for repeats and\
                         repeat element psuedogenomes.\
                         Example /data/annotation/mm9/setup''')
parser.add_argument('--nfragmentsfile1', action='store',
                    dest='nfragmentsfile1', metavar='nfragmentsfile1',
                    default='./repnames_nfragments.txt',
                    help='''File that saves the number\
                         of fragments processed per repname.
                         Default ./repnames_nfragments.txt''')
parser.add_argument('--gaplength', action='store', dest='gaplength',
                    metavar='gaplength', default='200', type=int,
                    help='Length of the spacer used to build\
                         repeat pseudogenomes.  Default 200')
parser.add_argument('--flankinglength', action='store', dest='flankinglength',
                    metavar='flankinglength', default='25', type=int,
                    help='Length of the flanking region adjacent to the repeat\
                         element that is used to build repeat pseudogenomes.\
                         The flanking length should be set according to the\
                         length of your reads.  Default 25')
args = parser.parse_args()

# parameters and paths specified in args_parse
gapl = args.gaplength
flankingl = args.flankinglength
annotation_file = args.annotation_file
genomefasta = args.genomefasta
setup_folder = args.setup_folder
nfragmentsfile1 = args.nfragmentsfile1

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
idxgenome = {}
lgenome = {}
genome = {}
allchrs = g.keys()
for k, chr in enumerate(allchrs):
    genome[chr] = str(g[chr].seq)
    lgenome[chr] = len(genome[chr])
    idxgenome[chr] = k

# Build a bedfile of repeatcoordinates to use by RepEnrich region_sorter
repeat_elements = []
fout = open(os.path.join(setup_folder, 'repnames.bed'), 'w')
rep_chr = {}
rep_start = {}
rep_end = {}
fin = import_text(annotation_file, ' ')
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
fin.close()
fout.close()

# sort repeat_elements and print them in repgenomes_key.txt
with open(os.path.join(setup_folder, 'repgenomes_key.txt'), 'w') as fout:
    for i, repeat in enumerate(sorted(repeat_elements)):
        fout.write('\t'.join([repeat, str(i)]) + '\n')

# generate spacer for pseudogenomes
spacer = ''.join(['N' for i in range(gapl)])

# save file with number of fragments processed per repname. Not used ????
with open(os.path.join(setup_folder, nfragmentsfile1), "w") as fout1:
    for repname in rep_chr:
        rep_chr_current = rep_chr[repname]
        fout1.write(str(len(rep_chr[repname])) + "\t" + repname + '\n')

# generate metagenomes and save them to FASTA files for bowtie build
nrepgenomes = len(rep_chr.keys())
for repname in rep_chr.keys():
    metagenome = ''
    rep_chr_current = rep_chr[repname]
    rep_start_current = rep_start[repname]
    rep_end_current = rep_end[repname]
    for i in range(len(rep_chr[repname])):
        try:
            chr = rep_chr_current[i]
            rstart = max(int(rep_start_current[i]) - flankingl, 0)
            rend = min(int(rep_end_current[i]) + flankingl, int(lgenome[chr])-1)
            metagenome = metagenome + spacer + genome[chr][rstart:(rend+1)]
        except KeyError:
            print("Unrecognised Chromosome: " + chr)
            pass
    # Convert metagenome to SeqRecord object (required by SeqIO.write)
    record = SeqRecord(Seq(metagenome), id='repname',
                       name='', description='')
    fastafilename = os.path.join(setup_folder, repname + '.fa')
    SeqIO.write(record, fastafilename, "fasta")
    command = shlex.split('bowtie-build -f ' + fastafilename + ' ' +
                          setup_folder + os.path.sep + repname)
    p = subprocess.Popen(command).communicate()
