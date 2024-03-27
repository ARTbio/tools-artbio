import argparse
import csv
import os
import shlex
import subprocess
import sys
from collections import defaultdict

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

# define a csv reader that reads space deliminated files
print('Preparing for analysis using RepEnrich...')
csv.field_size_limit(sys.maxsize)


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


def run_bowtie(metagenome, fastqfile, folder):
    output_file = os.path.join(folder, f"{metagenome}.bowtie")
    command = shlex.split(f"bowtie {b_opt} {metagenome} {fastqfile}")
    with open(output_file, 'w') as stdout:
        return subprocess.Popen(command, stdout=stdout)


# build dictionaries to convert repclass and rep families
repeatclass, repeatfamily = {}, {}
repeats = import_text(annotation_file, ' ')
for repeat in repeats:
    classfamily = []
    classfamily = repeat[10].split('/')
    matching_repeat = repeat[9].translate(str.maketrans('()/', '___'))
    repeatclass[matching_repeat] = classfamily[0]
    if len(classfamily) == 2:
        repeatfamily[matching_repeat] = classfamily[1]
    else:
        repeatfamily[matching_repeat] = classfamily[0]

# build list of repeats initializing dictionaries for downstream analysis
repgenome_path = 'repgenomes_key.txt'
reptotalcounts = {line.split('\t')[0]: 0 for line in open(repgenome_path)}
fractionalcounts = {line.split('\t')[0]: 0 for line in open(repgenome_path)}
classtotalcounts = {
    repeatclass[line.split('\t')[0]]: 0 for line in open(repgenome_path)
    if line.split('\t')[0] in repeatclass
}
classfractionalcounts = {
    repeatclass[line.split('\t')[0]]: 0 for line in open(repgenome_path)
    if line.split('\t')[0] in repeatclass
}
familytotalcounts = {
    repeatfamily[line.split('\t')[0]]: 0 for line in open(repgenome_path)
    if line.split('\t')[0] in repeatfamily
}
familyfractionalcounts = {
    repeatfamily[line.split('\t')[0]]: 0 for line in open(repgenome_path)
    if line.split('\t')[0] in repeatfamily
}
reptotalcounts_simple = {
    (simple_repeat if line.split('\t')[0] in repeatfamily
     and repeatfamily[line.split('\t')[0]] == simple_repeat
     else line.split('\t')[0]): 0 for line in open(repgenome_path)
}
repeat_key = {line.split('\t')[0]: int(line.split('\t')[1]) for line in open(
    repgenome_path)}
rev_repeat_key = {
    int(line.split('\t')[1]): line.split('\t')[0] for line in open(
        repgenome_path)}
repeat_list = [line.split('\t')[0] for line in open(repgenome_path)]

# map the repeats to the pseudogenomes
sorted_bowtie = 'sorted_bowtie'
os.makedirs(sorted_bowtie, exist_ok=True)

# unique mapper counting
cmd = f"bedtools bamtobed -i {unique_mapper_bam} | \
        bedtools coverage -b stdin -a  repnames.bed"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
bedtools_counts = p.communicate()[0].decode().rstrip('\r\n').split('\n')
counts = {}
sumofrepeatreads = 0
for line in bedtools_counts:
    line = line.split('\t')
    if not str(repeat_key[line[3]]) in counts:
        counts[str(repeat_key[line[3]])] = 0
    counts[str(repeat_key[line[3]])] += int(line[4])
    sumofrepeatreads += int(line[4])
print(f"Identified {sumofrepeatreads} unique reads that mapped to repeats.")

if not paired_end:
    folder_pair1 = 'pair1_bowtie'
    os.makedirs(folder_pair1, exist_ok=True)
    processes = []
    ticker = 0
    for metagenome in repeat_list:
        processes.append(run_bowtie(metagenome, fastqfile_1, folder_pair1))
        ticker += 1
        if ticker == cpus:
            for p in processes:
                p.communicate()
            ticker = 0
            processes = []
    for p in processes:
        p.communicate()
    for metagenome in repeat_list:
        file1 = os.path.join(folder_pair1, f"{metagenome}.bowtie")
        fileout = os.path.join(sorted_bowtie, f"{metagenome}.bowtie")
        collector = {}
        with open(fileout, 'w') as output:
            with open(file1) as input:
                for line in input:
                    collector[line.split("/")[0]] = 1
            for key in sorted(collector):
                output.write(f"{key}\n")
else:
    folder_pair1 = 'pair1_bowtie'
    folder_pair2 = 'pair2_bowtie'
    os.makedirs(folder_pair1, exist_ok=True)
    os.makedirs(folder_pair2, exist_ok=True)
    ps, psb, ticker = [], [], 0
    for metagenome in repeat_list:
        ps.append(run_bowtie(metagenome, fastqfile_1, folder_pair1))
        ticker += 1
        if fastqfile_2 != 'none':
            psb.append(run_bowtie(metagenome, fastqfile_2, folder_pair2))
            ticker += 1
        if ticker >= cpus:
            for p in ps:
                p.communicate()
            for p in psb:
                p.communicate()
            ticker = 0
            ps = []
            psb = []
    for p in ps:
        p.communicate()
    for p in psb:
        p.communicate()
    # collect read pair info
    for metagenome in repeat_list:
        file1 = folder_pair1 + os.path.sep + metagenome + '.bowtie'
        file2 = folder_pair2 + os.path.sep + metagenome + '.bowtie'
        fileout = sorted_bowtie + os.path.sep + metagenome + '.bowtie'
        collector = {}
        with open(fileout, 'w') as output:
            with open(file1) as input:
                for line in input:
                    collector[line.split("/")[0]] = 1
            with open(file2) as input:
                for line in input:
                    collector[line.split("/")[0]] = 1
            for key in sorted(collector):
                output.write(f"{key}\n")

# build a file of repeat keys for all reads
sumofrepeatreads = 0
readid = defaultdict(dict)

for rep in repeat_list:
    for line in open( f"{os.path.join(sorted_bowtie, rep)}.bowtie"):
        readid[line.split()[0]] += f"{repeat_key[rep]},"

for subfamilies in readid.values():
    if subfamilies not in counts:
        counts[subfamilies] = 0
    counts[subfamilies] += 1
    sumofrepeatreads += 1

print(f'Identified {sumofrepeatreads} reads that mapped to \
        repeats for unique and multimappers.')
print("Conducting final calculations...")

# building the total counts for repeat element enrichment...
for x in counts.keys():
    count = counts[x]
    x = x.strip(',').split(',')
    for i in x:
        reptotalcounts[rev_repeat_key[int(i)]] += int(count)
# building the fractional counts for repeat element enrichment...
for x in counts.keys():
    count = counts[x]
    x = x.strip(',').split(',')
    splits = len(x)
    for i in x:
        fractionalcounts[rev_repeat_key[int(i)]] += float(
            numpy.divide(float(count), float(splits)))
# building categorized table of repeat element enrichment...
repcounts = {}
repcounts['other'] = 0
for key in counts.keys():
    key_list = key.strip(',').split(',')
    repname = ''
    for i in key_list:
        repname = os.path.join(repname, rev_repeat_key[int(i)])
    repcounts[repname] = counts[key]
# building the total counts for class enrichment...
for key in reptotalcounts.keys():
    classtotalcounts[repeatclass[key]] += reptotalcounts[key]
# building total counts for family enrichment...
for key in reptotalcounts.keys():
    familytotalcounts[repeatfamily[key]] += reptotalcounts[key]
# building unique counts table
repcounts2 = {}
for rep in repeat_list:
    if "/" + rep in repcounts:
        repcounts2[rep] = repcounts["/" + rep]
    else:
        repcounts2[rep] = 0
# building the fractionalcounts counts for class enrichment...
for key in fractionalcounts.keys():
    classfractionalcounts[repeatclass[key]] += fractionalcounts[key]
# building fractional counts for family enrichment...
for key in fractionalcounts.keys():
    familyfractionalcounts[repeatfamily[key]] += fractionalcounts[key]

# print output to file of the categorized counts and total overlapping counts:
print('Writing final output...')
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
