import argparse

from collections import defaultdict

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('-bam', '--bam', dest='bams', required=True,
                            nargs='+', help='input BAM file')
    the_parser.add_argument('-o', '--output', dest='distribution', required=True,
                            help='tabular output for mapq distribution')
    args = the_parser.parse_args()
    return args


def collect_mapq(bam, out):
    samfile = pysam.AlignmentFile(bam, "rb")
    mapq_dict = defaultdict(int)
    for read in samfile:
        mapq_dict[read.mapping_quality] += 1
    with open(out, 'w') as out:
        out.write('mapq\tnumber_of_alignments')
        for quality in sorted(mapq_dict):
            out.write(f"{quality}\t{mapq_dict[quality]}")
    return mapq_dict


def main(bam, out):
    results_dict = collect_mapq(bam, out)


if __name__ == "__main__":
    args = Parser()
    main(args.bam, args.out)
