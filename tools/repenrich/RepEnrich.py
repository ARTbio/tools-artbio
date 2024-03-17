import argparse
import csv
import os
import shlex
import subprocess
import sys

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
parser.add_argument('--outputfolder', action='store', metavar='outputfolder',
                    help='Folder that will contain results. Should be the\
                          same as the one used for RepEnrich_setup.\
                          Example: ./outputfolder')
parser.add_argument('--outputprefix', action='store', metavar='outputprefix',
                    help='Prefix name for Repenrich output files.')
parser.add_argument('--setup_folder', action='store', metavar='setup_folder',
                    help='Folder produced by RepEnrich_setup which contains\
                    repeat element pseudogenomes.')
parser.add_argument('--fastqfile', action='store', metavar='fastqfile',
                    help='File of fastq reads mapping to multiple\
                          locations. Example: /data/multimap.fastq')
parser.add_argument('--alignment_bam', action='store', metavar='alignment_bam',
                    help='Bam alignments of unique mapper reads.')
parser.add_argument('--pairedend', action='store', dest='pairedend',
                    default='FALSE',
                    help='Change to TRUE for paired-end fastq files.\
                          Default FALSE')
parser.add_argument('--collapserepeat', action='store', dest='collapserepeat',
                    metavar='collapserepeat', default='Simple_repeat',
                    help='Use this option to generate a collapsed repeat\
                          type. Uncollapsed output is generated in addition to\
                          collapsed repeat type. Simple_repeat is default to\
                          simplify downstream analysis. You can change the\
                          default to another repeat name to collapse a\
                          separate specific repeat instead or if the name of\
                          Simple_repeat is different for your organism.\
                          Default: "Simple_repeat"')
parser.add_argument('--fastqfile2', action='store', dest='fastqfile2',
                    metavar='fastqfile2', default='none',
                    help='fastqfile #2 when using paired-end option.\
                          Default none')
parser.add_argument('--cpus', action='store', dest='cpus', metavar='cpus',
                    default="1", type=int,
                    help='Number of CPUs. The more cpus the\
                          faster RepEnrich performs. Default: "1"')
parser.add_argument('--allcountmethod', action='store', dest='allcountmethod',
                    metavar='allcountmethod', default="FALSE",
                    help='By default the script only outputs the fraction\
                          counts.\
                          Other options\
                          include the unique count method, a conservative\
                          count, and the total count method, a liberal\
                          counting strategy.\
                          Our evaluation of simulated data indicated\
                          that fraction counting is the best method.\
                          Default = FALSE, change to TRUE')

args = parser.parse_args()

# parameters
annotation_file = args.annotation_file
outputfolder = args.outputfolder
outputfile_prefix = args.outputprefix
setup_folder = args.setup_folder
repeat_bed = os.path.join(setup_folder, 'repnames.bed')
unique_mapper_bam = args.alignment_bam
fastqfile_1 = args.fastqfile
fastqfile_2 = args.fastqfile2
cpus = args.cpus
b_opt = "-k1 -p 1 --quiet"
simple_repeat = args.collapserepeat
paired_end = args.pairedend
allcountmethod = args.allcountmethod


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


def import_text(filename, separator):
    for line in csv.reader(open(filename), delimiter=separator,
                           skipinitialspace=True):
        if line:
            yield line


# build dictionaries to convert repclass and rep families
repeatclass, repeatfamily = {}, {}
repeats = import_text(annotation_file, ' ')
# skip three first lines of the iterator
for line in range(3):
    next(repeats)
for repeat in repeats:
    classfamily = []
    classfamily = repeat[10].split('/')
    matching_repeat = repeat[9].translate(str.maketrans('()/', '___'))
    repeatclass[matching_repeat] = classfamily[0]
    if len(classfamily) == 2:
        repeatfamily[matching_repeat] = classfamily[1]
    else:
        repeatfamily[matching_repeat] = classfamily[0]


# build list of repeats initializing dictionaries for downstream analysis'
repgenome_path = os.path.join(setup_folder, 'repgenomes_key.txt')
reptotalcounts = {line[0]: 0 for line in import_text(repgenome_path, '\t')}
fractionalcounts = {line[0]: 0 for line in import_text(repgenome_path, '\t')}
classtotalcounts = {repeatclass[line[0]]: 0 for line in import_text(
    repgenome_path, '\t') if line[0] in repeatclass}
classfractionalcounts = {repeatclass[line[0]]: 0 for line in import_text(
    repgenome_path, '\t') if line[0] in repeatclass}
familytotalcounts = {repeatfamily[line[0]]: 0 for line in import_text(
    repgenome_path, '\t') if line[0] in repeatfamily}
familyfractionalcounts = {repeatfamily[line[0]]: 0 for line in import_text(
    repgenome_path, '\t') if line[0] in repeatfamily}
reptotalcounts_simple = {(simple_repeat if line[0] in repeatfamily and
                          repeatfamily[line[0]] == simple_repeat else
                          line[0]): 0 for line in import_text(
                              repgenome_path, '\t')}
repeat_key = {line[0]: int(line[1]) for line in import_text(
    repgenome_path, '\t')}
rev_repeat_key = {int(line[1]): line[0] for line in import_text(
    repgenome_path, '\t')}
repeat_list = [line[0] for line in import_text(repgenome_path, '\t')]

# map the repeats to the psuedogenomes:
if not os.path.exists(outputfolder):
    os.mkdir(outputfolder)

# Conduct the regions sorting
fileout = os.path.join(outputfolder, f"{outputfile_prefix}_regionsorter.txt")
command = shlex.split(f"coverageBed -abam {unique_mapper_bam} -b \
                        {os.path.join(setup_folder, 'repnames.bed')}")
with open(fileout, 'w') as stdout:
    subprocess.run(command, stdout=stdout, check=True)

counts = {}
sumofrepeatreads = 0
with open(fileout) as filein:
    for line in filein:
        line = line.split('\t')
        if not str(repeat_key[line[3]]) in counts:
            counts[str(repeat_key[line[3]])] = 0
        counts[str(repeat_key[line[3]])] += int(line[4])
        sumofrepeatreads += int(line[4])
    print(f"Identified {sumofrepeatreads} unique reads that \
            mapped to repeats.")


def run_bowtie(metagenome, fastqfile, folder):
    metagenomepath = os.path.join(setup_folder, metagenome)
    output_file = os.path.join(folder, f"{metagenome}.bowtie")
    command = shlex.split(f"bowtie {b_opt} {metagenomepath} {fastqfile}")
    with open(output_file, 'w') as stdout:
        return subprocess.Popen(command, stdout=stdout)


if paired_end == 'FALSE':
    folder_pair1 = os.path.join(outputfolder, 'pair1_bowtie')
    os.makedirs(folder_pair1, exist_ok=True)

    print("Processing repeat pseudogenomes...")
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

    # Combine the output from both read pairs:
    print('Sorting and combining the output for both read pairs....')
    sorted_bowtie = os.path.join(outputfolder, 'sorted_bowtie')
    os.makedirs(sorted_bowtie, exist_ok=True)
    for metagenome in repeat_list:
        file1 = os.path.join(folder_pair1, f"{metagenome}.bowtie")
        fileout = os.path.join(sorted_bowtie, f"{metagenome}.bowtie")
        with open(fileout, 'w') as stdout:
            p1 = subprocess.Popen(['cat', file1], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['cut', '-f1'], stdin=p1.stdout,
                                  stdout=subprocess.PIPE)
            p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin=p2.stdout,
                                  stdout=subprocess.PIPE)
            p4 = subprocess.Popen(['sort'], stdin=p3.stdout,
                                  stdout=subprocess.PIPE)
            p5 = subprocess.Popen(['uniq'], stdin=p4.stdout, stdout=stdout)
            p5.communicate()
        stdout.close()
    print('completed ...')
else:
    folder_pair1 = os.path.join(outputfolder, 'pair1_bowtie')
    folder_pair2 = os.path.join(outputfolder, 'pair2_bowtie')
    os.makedirs(folder_pair1, exist_ok=True)
    os.makedirs(folder_pair2, exist_ok=True)

    print("Processing repeat pseudogenomes...")
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

    # combine the output from both read pairs:
    print('Sorting and combining the output for both read pairs...')
    if not os.path.exists(outputfolder + os.path.sep + 'sorted_bowtie'):
        os.mkdir(outputfolder + os.path.sep + 'sorted_bowtie')
    sorted_bowtie = outputfolder + os.path.sep + 'sorted_bowtie'
    for metagenome in repeat_list:
        file1 = folder_pair1 + os.path.sep + metagenome + '.bowtie'
        file2 = folder_pair2 + os.path.sep + metagenome + '.bowtie'
        fileout = sorted_bowtie + os.path.sep + metagenome + '.bowtie'
        with open(fileout, 'w') as stdout:
            p1 = subprocess.Popen(['cat', file1, file2],
                                  stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['cut', '-f1', "-d "], stdin=p1.stdout,
                                  stdout=subprocess.PIPE)
            p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin=p2.stdout,
                                  stdout=subprocess.PIPE)
            p4 = subprocess.Popen(['sort'], stdin=p3.stdout,
                                  stdout=subprocess.PIPE)
            p5 = subprocess.Popen(['uniq'], stdin=p4.stdout, stdout=stdout)
            p5.communicate()
        stdout.close()
    print('completed ...')

# build a file of repeat keys for all reads
print('Writing and processing intermediate files...')
sorted_bowtie = os.path.join(outputfolder, 'sorted_bowtie')
sumofrepeatreads = 0
readid = {}

for rep in repeat_list:
    for data in import_text(
            f"{os.path.join(sorted_bowtie, rep)}.bowtie", '\t'):
        readid[data[0]] = ''

for rep in repeat_list:
    for data in import_text(
            f"{os.path.join(sorted_bowtie, rep)}.bowtie", '\t'):
        readid[data[0]] += f"{repeat_key[rep]},"

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
    x = x.strip(',')    .split(',')
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

print('Writing final output and removing intermediate files...')
# print output to file of the categorized counts and total overlapping counts:
if allcountmethod == "TRUE":
    fout1 = open(outputfolder + os.path.sep + outputfile_prefix
                 + '_total_counts.txt', 'w')
    for key in sorted(reptotalcounts.keys()):
        fout1.write(str(key) + '\t' + repeatclass[key] + '\t' +
                    repeatfamily[key] + '\t' + str(reptotalcounts[key])
                    + '\n')
    fout2 = open(outputfolder + os.path.sep + outputfile_prefix
                 + '_class_total_counts.txt', 'w')
    for key in sorted(classtotalcounts.keys()):
        fout2.write(str(key) + '\t' + str(classtotalcounts[key]) + '\n')
    fout3 = open(outputfolder + os.path.sep + outputfile_prefix
                 + '_family_total_counts.txt', 'w')
    for key in sorted(familytotalcounts.keys()):
        fout3.write(str(key) + '\t' + str(familytotalcounts[key]) + '\n')
    fout4 = open(outputfolder + os.path.sep + outputfile_prefix +
                 '_unique_counts.txt', 'w')
    for key in sorted(repcounts2.keys()):
        fout4.write(str(key) + '\t' + repeatclass[key] + '\t' +
                    repeatfamily[key] + '\t' + str(repcounts2[key]) + '\n')
        fout5 = open(outputfolder + os.path.sep + outputfile_prefix
                     + '_class_fraction_counts.txt', 'w')
    for key in sorted(classfractionalcounts.keys()):
        fout5.write(str(key) + '\t' + str(classfractionalcounts[key]) + '\n')
    fout6 = open(outputfolder + os.path.sep + outputfile_prefix +
                 '_family_fraction_counts.txt', 'w')
    for key in sorted(familyfractionalcounts.keys()):
        fout6.write(str(key) + '\t' + str(familyfractionalcounts[key]) + '\n')
    fout7 = open(outputfolder + os.path.sep + outputfile_prefix
                 + '_fraction_counts.txt', 'w')
    for key in sorted(fractionalcounts.keys()):
        fout7.write(str(key) + '\t' + repeatclass[key] + '\t' +
                    repeatfamily[key] + '\t' + str(int(fractionalcounts[key]))
                    + '\n')
        fout1.close()
    fout2.close()
    fout3.close()
    fout4.close()
    fout5.close()
    fout6.close()
    fout7.close()
else:
    fout1 = open(outputfolder + os.path.sep + outputfile_prefix +
                 '_class_fraction_counts.txt', 'w')
    for key in sorted(classfractionalcounts.keys()):
        fout1.write(str(key) + '\t' + str(classfractionalcounts[key]) + '\n')
    fout2 = open(outputfolder + os.path.sep + outputfile_prefix +
                 '_family_fraction_counts.txt', 'w')
    for key in sorted(familyfractionalcounts.keys()):
        fout2.write(str(key) + '\t' + str(familyfractionalcounts[key]) + '\n')
    fout3 = open(outputfolder + os.path.sep + outputfile_prefix +
                 '_fraction_counts.txt', 'w')
    for key in sorted(fractionalcounts.keys()):
        fout3.write(str(key) + '\t' + repeatclass[key] + '\t' +
                    repeatfamily[key] + '\t' + str(int(fractionalcounts[key]))
                    + '\n')
    fout1.close()
    fout2.close()
    fout3.close()
