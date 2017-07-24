#!/usr/bin/env python
import argparse

import pysam


def Parser():
    parser = argparse.ArgumentParser(description='miRNAs counts and coverages')
    parser.add_argument('-a', '--alignment', metavar='FILE', type=str,
                        dest='alignment_file', help='Alignment bam file')
    parser.add_argument('--gff', metavar='FILE', type=str, dest='gff_file',
                        help='GFF3 describing both pre-miRNAs\
                              and mature miRNAs')
    parser.add_argument('-q', '--quality_threshold', type=int,
                        dest='quality_threshold',
                        help='Quality threshold for coverage (default=10)',
                        default=10)
    parser.add_argument('-p', '--pre_mirs', type=str, dest='pre_mirs',
                        help='pre-miRNAs count file path', metavar='FILE')
    parser.add_argument('-m', '--mirs', type=str, dest='mirs',
                        help='mature miRNA count file path', metavar='FILE')
    parser.add_argument('--lattice', metavar='FILE', type=str, dest='lattice',
                        help='Output file for the lattice dataframe.')
    args = parser.parse_args()
    return args


def get_pre_mir_counts(bamfile):
    """
    Takes a AlignmentFile object and returns a dictionary of counts for reads
    aligning with pre_mirs (as keys)
    """
    count = dict()
    for ref_name in bamfile.references:
        count[ref_name] = bamfile.count(reference=ref_name)
    return count


def get_pre_mir_coverage(bamfile, quality=10):
    """
    Takes a AlignmentFile object and returns a dictionary of lists
    of coverage along the coordinates of pre_mirs (as keys)
    """
    coverage = dict()
    for ref_name, ref_len in zip(bamfile.references, bamfile.lengths):
        coverage[ref_name] = bamfile.count_coverage(reference=ref_name,
                                                    start=0, end=ref_len,
                                                    quality_threshold=quality)
        """ Add the 4 coverage values """
        coverage[ref_name] = [sum(x) for x in
                              zip(*coverage[ref_name])]
    return coverage


def get_mir_counts(bamfile, gff_file):
    """
    Takes a AlignmentFile and a gff file and computes for
    each 'miRNA' region of the gff the number of reads that hit it
    returns a dict[mir_name] = count
    """
    counts = dict()
    for line in open(gff_file, 'r'):
        if line[0] != '#':
            gff_fields = line[:-1].split("\t")
            if gff_fields[2] == 'miRNA':
                mir_name = gff_fields[0]
                premir_name = gff_fields[8].split('=')[-1]
                mir_start = int(gff_fields[3])
                mir_end = int(gff_fields[4])
                # GFF is 1-based, pysam is 0-based.
                counts[mir_name] = bamfile.count(reference=premir_name,
                                                 start=mir_start-1,
                                                 end=mir_end-1)
    return counts


def write_dataframe_coverage(countdict, outfile):
    """
    Takes a dict[pre_mir reference name] = [coverage list]
    and writes a dataframe with columns:
    <gene_type name>, offset, normoffset, counts and normcounts
    in the outfile
    """
    F = open(outfile, 'w')
    F.write('Mir_hairpin\tOffset\tNorm_offset\tCount\tNorm_count\n')
    for ref in sorted(countdict):
        """
        For each reference name in mirs,
        write the coverage of each of its positions
        """
        maximum = max(countdict[ref])
        reference_length = len(countdict[ref])
        for pos, c in enumerate(countdict[ref]):
            """ Compute and write value for each reference position"""
            F.write('%s\t%s\t%s\t%s\t%s\n' % (ref, str(pos + 1),
                    str(float(pos+1)/reference_length), str(float(c)),
                    str(float(c)/maximum) if maximum != 0 else '0'))
    F.close()


def write_counts(countdict, outfile):
    """
    Takes a dict[<gene_type name>]=count and
    writes a count table
    """
    F = open(outfile, 'w')
    F.write('Gene\tCounts\n')
    for gene in sorted(countdict):
        F.write('%s\t%s\n' % (gene, str(countdict[gene])))
    F.close()


def main():
    args = Parser()
    bamfile = pysam.AlignmentFile(args.alignment_file, 'rb', check_sq=False)
    if args.pre_mirs:
        pre_mirs = get_pre_mir_counts(bamfile)
        write_counts(pre_mirs, args.pre_mirs)
        if args.lattice:
            pre_mirs_coverage = get_pre_mir_coverage(bamfile,
                                                     args.quality_threshold)
            write_dataframe_coverage(pre_mirs_coverage, args.lattice)
    if args.mirs:
        mirs = get_mir_counts(bamfile, args.gff_file)
        write_counts(mirs, args.mirs)


if __name__ == '__main__':
    main()
