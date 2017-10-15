#!/usr/bin/env python
import argparse
from collections import defaultdict

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--output', nargs='+', action='store', type=str,
                            help='Count tables')
    the_parser.add_argument('--polarity',
                            choices=["sense", "antisense", "both"],
                            help="forward, reverse or both forward an\
                                reverse reads are counted")
    the_parser.add_argument('--alignments', nargs='+',
                            help="bam alignments files")
    the_parser.add_argument('--labels', nargs='+', help="Alignments labels")
    the_parser.add_argument('--number',
                            choices=["unique", "multiple"],
                            help="output is a single table or multiple tables")
    args = the_parser.parse_args()
    return args


def get_counts(bamfile, polarity="both"):
    """
    Takes an AlignmentFile object and returns a dictionary of counts for sense,
    antisense, or both sense and antisense reads aligning to the bam references
    """
    def filter_sense_read(read):
        if read.is_reverse:
            return 0
        else:
            return 1

    def filter_antisense_read(read):
        if read.is_reverse:
            return 1
        else:
            return 0

    counts = defaultdict(int)
    for ref_name in bamfile.references:
        counts[ref_name] = 0
    if polarity == "both":
        for ref_name in bamfile.references:
            counts[ref_name] = bamfile.count(reference=ref_name)
    if polarity == "sense":
        for ref_name in bamfile.references:
            for read in bamfile.fetch(ref_name):
                counts[ref_name] += filter_sense_read(read)
    if polarity == "antisense":
        for ref_name in bamfile.references:
            for read in bamfile.fetch(ref_name):
                counts[ref_name] += filter_antisense_read(read)
    return counts


def writetable(diclist, labels, output, number):
    ''' diclist is a list of count dictionnaries '''
    countlists = []
    for dic in diclist:
        counts = sorted(dic.items())
        counts = [j for (i, j) in counts]
        countlists.append(counts)
    if number == "unique":
        out = open("outputdir/table.tabular", "w")
        out.write("gene\t%s\n" % "\t".join(labels))
        for countline in zip(sorted(diclist[0]), *countlists):
            line = [str(i) for i in countline]
            out.write("%s\n" % "\t".join(line))
        out.close()
    else:
        for i, (dic, label) in enumerate(zip(diclist, labels)):
            out = open("outputdir/table" + str(i) + ".tabular", "w")
            out.write("gene\t%s\n" % label)
            for gene in sorted(dic):
                out.write("%s\t%s\n" % (gene, dic[gene]))
            out.close()


def main(alignments, labels, polarity, output, number):
    diclist = []
    for file in alignments:
        bam_object = pysam.AlignmentFile(file, 'rb')
        diclist.append(get_counts(bam_object, polarity=polarity))
    writetable(diclist, labels, output, number)


if __name__ == "__main__":
    args = Parser()
    main(args.alignments, args.labels, args.polarity, args.output, args.number)
