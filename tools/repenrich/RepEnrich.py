import argparse
import csv
import os
import shlex
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy


parser = argparse.ArgumentParser(description='''
             Repenrich aligns reads to Repeat Elements pseudogenomes\
             and counts aligned reads. RepEnrich_setup must be run\
             before its use''')
parser.add_argument('--annotation_file', action='store',
                    metavar='annotation_file',
                    help='RepeatMasker.org annotation file for your\
                          organism. The file may be downloaded from\
                          RepeatMasker.org. E.g. hg19_repeatmasker.txt')
parser.add_argument('--alignment_bam', action='store', metavar='alignment_bam',
                    help='Bam alignments of unique mapper reads.')
parser.add_argument('--fastqfile', action='store', metavar='fastqfile',
                    help='File of fastq reads mapping to multiple\
                          locations. Example: /data/multimap.fastq')
parser.add_argument('--fastqfile2', action='store', dest='fastqfile2',
                    metavar='fastqfile2', default='',
                    help='fastqfile #2 when using paired-end option.\
                          Default none')
parser.add_argument('--cpus', action='store', dest='cpus', metavar='cpus',
                    default="1", type=int,
                    help='Number of CPUs. The more cpus the\
                          faster RepEnrich performs. Default: "1"')

args = parser.parse_args()

# parameters
annotation_file = args.annotation_file
repeat_bed = 'repnames.bed'
unique_mapper_bam = args.alignment_bam
fastqfile_1 = args.fastqfile
fastqfile_2 = args.fastqfile2
cpus = args.cpus
b_opt = "-k1 -p 1 --quiet"
# Change if simple repeats are differently annotated in your organism
simple_repeat = "Simple_repeat"
if args.fastqfile2:
    paired_end = True
else:
    paired_end = False

# check that the programs we need are available
try:
    subprocess.call(shlex.split("coverageBed -h"),
                    stdout=open(os.devnull, 'wb'),
                    stderr=open(os.devnull, 'wb'))
    subprocess.call(shlex.split("bowtie --version"),
                    stdout=open(os.devnull, 'wb'),
                    stderr=open(os.devnull, 'wb'))
except OSError:
    print("Error: Bowtie or bedtools not loaded")
    raise


def starts_with_numerical(list):
    try:
        if len(list) == 0:
            return False
        int(list[0])
        return True
    except ValueError:
        return False


# define a text importer for .out/.txt format of repbase
def import_text(filename, separator):
    csv.field_size_limit(sys.maxsize)
    file = csv.reader(open(filename), delimiter=separator,
                      skipinitialspace=True)
    return [line for line in file if starts_with_numerical(line)]


def run_bowtie(args):
    metagenome, fastqfile, folder = args
    output_file = os.path.join(folder, f"{metagenome}.bowtie")
    command = shlex.split(f"bowtie {b_opt} {metagenome} {fastqfile}")
    with open(output_file, 'w') as stdout:
        return subprocess.Popen(command, stdout=stdout)


repeat_key = {line.split('\t')[0]: int(line.split('\t')[1]) for line in open(
    'repeatIDs.txt')}
rev_repeat_key = {v: k for k, v in repeat_key.items()}
repeat_list = [line.split('\t')[0] for line in open('repeatIDs.txt')]

# map the repeats to the pseudogenomes
sorted_bowtie = 'sorted_bowtie'
os.makedirs(sorted_bowtie, exist_ok=True)

# unique mapper counting
cmd = f"bedtools bamtobed -i {unique_mapper_bam} | \
        bedtools coverage -b stdin -a  repnames.bed"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
bedtools_counts = p.communicate()[0].decode().rstrip('\r\n').split('\n')

# parse bedtools output
counts = defaultdict(int)
sumofrepeatreads = 0
for line in bedtools_counts:
    line = line.split('\t')
    counts[str(repeat_key[line[3]])] += int(line[4])
    sumofrepeatreads += int(line[4])
print(f"Identified {sumofrepeatreads} unique reads that mapped to repeats.")

# Analyse multimappers
if not paired_end:
    folder_pair1 = 'pair1_bowtie'
    os.makedirs(folder_pair1, exist_ok=True)
    args_list = [(metagenome, fastqfile_1, folder_pair1) for
                 metagenome in repeat_list]

    with ProcessPoolExecutor(max_workers=cpus) as executor:
        executor.map(run_bowtie, args_list)
    # Consolidate the output files
    for metagenome in repeat_list:
        file1 = os.path.join(folder_pair1, f"{metagenome}.bowtie")
        fileout = os.path.join(sorted_bowtie, f"{metagenome}.bowtie")
        collector = {}
        with open(file1) as input_file:
            for line in input_file:
                collector[line.split("/")[0]] = 1
        with open(fileout, 'w') as output:
            for key in sorted(collector):
                output.write(f"{key}\n")
else:
    folder_pair1 = 'pair1_bowtie'
    folder_pair2 = 'pair2_bowtie'
    os.makedirs(folder_pair1, exist_ok=True)
    os.makedirs(folder_pair2, exist_ok=True)
    args_list = [(metagenome, fastqfile_1, folder_pair1) for
                 metagenome in repeat_list]
    args_list.extend([(metagenome, fastqfile_2, folder_pair2) for
                     metagenome in repeat_list])
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        executor.map(run_bowtie, args_list)
    # collect read pair info
    for metagenome in repeat_list:
        file1 = os.path.join(folder_pair1, f"{metagenome}.bowtie")
        file2 = os.path.join(folder_pair2, f"{metagenome}.bowtie")
        fileout = os.path.join(sorted_bowtie, f"{metagenome}.bowtie")
        collector = {}
        with open(file1) as input:
            for line in input:
                collector[line.split("/")[0]] = 1
        with open(file2) as input:
            for line in input:
                collector[line.split("/")[0]] = 1
        with open(fileout, 'w') as output:
            for key in sorted(collector):
                output.write(f"{key}\n")

# build dictionaries to convert repclass and rep families
repeatclass, repeatfamily = {}, {}
repeats = import_text(annotation_file, ' ')
for repeat in repeats:
    classfamily = repeat[10].split('/')
    matching_repeat = repeat[9].translate(str.maketrans('()/', '___'))
    repeatclass[matching_repeat] = classfamily[0]
    if len(classfamily) == 2:
        repeatfamily[matching_repeat] = classfamily[1]
    else:
        repeatfamily[matching_repeat] = classfamily[0]

# build list of repeats initializing dictionaries for downstream analysis
readid = defaultdict(list)
# counts dictionary already implemented above
reptotalcounts = defaultdict(int)
familytotalcounts = defaultdict(int)
classtotalcounts = defaultdict(int)
fractionalcounts = defaultdict(float)
familyfractionalcounts = defaultdict(float)
classfractionalcounts = defaultdict(float)
repcounts = {}
repcounts2 = {}
sumofrepeatreads = 0

# Loop through bowtie files and populate readid
for rep in repeat_list:
    with open(os.path.join(sorted_bowtie, f"{rep}.bowtie")) as file:
        for line in file:
            read, _ = line.split(maxsplit=1)  # Split only once
            readid[read].append(repeat_key[rep])

# Populate counts and sumofrepeatreads
for subfamilies in readid.values():
    subfamilies_str = ','.join(subfamilies)
    counts[subfamilies_str] += 1
    sumofrepeatreads += 1

print(f'Identified {sumofrepeatreads} reads that mapped to \
        repeats for unique and multimappers.')

# Populate reptotalcounts and fractionalcounts
for key, count in counts.items():
    key_list = key.split(',')
    for i in key_list:
        reptotalcounts[rev_repeat_key[int(i)]] += count
    splits = len(key_list)
    for i in key_list:
        fractionalcounts[rev_repeat_key[int(i)]] += count / splits

# Populate repcounts
for key, count in counts.items():
    key_list = key.split(',')
    repname = '/'.join(rev_repeat_key[int(i)] for i in key_list)
    repcounts[repname] = count

# Populate classtotalcounts and familytotalcounts
for key, value in reptotalcounts.items():
    classtotalcounts[repeatclass[key]] += value
    familytotalcounts[repeatfamily[key]] += value

# Populate repcounts2
for rep in repeat_list:
    repcounts2[rep] = repcounts.get('/' + rep, 0)

# Populate classfractionalcounts and familyfractionalcounts
for key, value in fractionalcounts.items():
    classfractionalcounts[repeatclass[key]] += value
    familyfractionalcounts[repeatfamily[key]] += value

# print output to files of the categorized counts and total overlapping counts:
with open("class_fraction_counts.tsv", 'w') as fout:
    for key in sorted(classfractionalcounts.keys()):
        fout.write(f"{key}\t{classfractionalcounts[key]}\n")

with open("family_fraction_counts.tsv", 'w') as fout:
    for key in sorted(familyfractionalcounts.keys()):
        fout.write(f"{key}\t{familyfractionalcounts[key]}\n")

with open("fraction_counts.tsv", 'w') as fout:
    for key in sorted(fractionalcounts.keys()):
        fout.write(f"{key}\t{repeatclass[key]}\t{repeatfamily[key]}\t"
                   f"{int(fractionalcounts[key])}\n")
