import argparse

import numpy

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('-bams', '--bams', dest='bams', required=True,
                            nargs='+', help='list of input BAM files')
    the_parser.add_argument('-bed', '--bed', dest='bed', required=True,
                            help='Coordinates of probes in a bed file')
    args = the_parser.parse_args()
    return args


def compute_coverage(bam, bed, quality=10):
    bam_object = pysam.AlignmentFile(bam, 'rb')
    bed_object = open(bed, 'r')
    coverage_column = []
    for line in bed_object:
        if line[0] == '#':
            continue
        fields = line[:-1].split('\t')
        chr = fields[0]
        start = fields[1]
        end = fields[2]
        coverage = bam_object.count_coverage(reference=chr,
                                             start=int(start)-1,
                                             stop=int(end),
                                             quality_threshold=quality)
        """ Add the 4 coverage values """
        coverage = [sum(x) for x in zip(*coverage)]
        coverage_column.append(numpy.mean(coverage))
    bed_object.close()
    return coverage_column


def main(bams, bed):
    column_dict = {}
    for i, bam in enumerate(bams):
        column_dict[i] = compute_coverage(bam, bed)
    F = open(bed, 'r')
    line_counter = 0
    for line in F:
        if line[0] == '#':
            continue
        prefix = line[:-1]
        crossline = []
        for col in sorted(column_dict):
            crossline.append(str(column_dict[col][line_counter]))
        line_counter += 1
        suffix = '\t'.join(crossline)
        print('%s\t%s' % (prefix, suffix))


if __name__ == "__main__":
    args = Parser()
    main(args.bams, args.bed)
