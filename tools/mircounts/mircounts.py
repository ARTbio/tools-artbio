#!/usr/bin/env python

import pysam
import sys
import optparse

def option_parsing():
    parser = optparse.OptionParser(description='Count miRNAs')
    parser.add_option('-p', '--pre_mirs', action='store_true', dest='pre_mirs', help='Count pre-miRNAs', default=False)
    parser.add_option('-m', '--mirs', action='store_true', dest='mirs', help='Count mature miRNAs', default=False)
    inputs = optparse.OptionGroup(parser, 'Inputs')
    inputs.add_option('--alignment', type='string', dest='alignment_file', help='Alignment tabular or sam file')
    inputs.add_option('--gff', type='string', dest='gff_file', help=' GFF3 describing both pre-mirs and mature mirs (to be discussed)', metavar='FILE')
    parser.add_option_group(inputs)
    outputs = optparse.OptionGroup(parser, 'Outputs')
    outputs.add_option('--pre_mirs_output', type='string', dest='output_pre_mirs', help='GFF3 describing the mature miRs', metavar='FILE', default='output_pre_mirs_count.tab')
    outputs.add_option('--mature_mirs_output', type='string', dest='output_mature_mirs', help='Reference genome', metavar='FILE', default='output_mirs_count.tab')
    outputs.add_option('--lattice', type='string', dest='lattice')
    parser.add_option_group(outputs)
    (options, args) = parser.parse_args()
    """ Check if the options were correctly passed """
    if len(args) > 0:
        parser.error('Wrong number of arguments')
    if (not options.alignment_file or not options.gff_file):
        parser.error("Missing file. Both '--alignment' and '--gff' files are needed")
    if (options.lattice and not options.mirs):
        parse.error("Can't output lattice dataframe if the option '-m' is not set")
    return options

def get_count(bamfile):
    """
    Takes a AlignmentFile object and and returns a count of reads
    hiting 'miRNA_primary_transcript' positions
    """
    count = dict()
    for reference in bamfile.header['SQ']:
        count[reference['SN']] = [bamfile.count(reference['SN'])]
    return count

def get_coverage(bamfile, gff_file):
    """
    Takes a AlignmentFile and a gff file and computes for
    each 'miRNA' region of the gff the number of reads that hit it
    and the coverage of each of its positions
    and returns a dictionnary : dict[reference name] = [polarity, coverage as [[],[],[],[]]]
    """
    counts = dict()
    try:
        gff = open(gff_file, 'r')
        for line in gff.readlines():
            if line[0] != '#':
                gff_fields = line[:-1].split("\t")
                if gff_fields[2] == 'miRNA':
                    counts[gff_fields[0]] = [bamfile.count(reference=gff_fields[0], start=int(gff_fields[3]),
                                                           end=int(gff_fields[4])),
                                             gff_fields[6],
                                             bamfile.count_coverage(reference=gff_fields[0], start=int(gff_fields[3]),
                                                                    end=int(gff_fields[4]))]
        gff.close()
    except IOError as e:
        sys.stderr.write("Error while reading file %s\n" % gff_file)
        sys.stderr.write("I/O error(%s): %s\n" % (e.errno, e.strerror))
    return counts

def write_dataframe(mirs, outfile, sample):
    """
    Takes a dictionnary dict[reference name] = [polarity, count,
                                                [[Hits of each position for 'A'],
                                                 [Hits of each position for 'C'],
                                                 [Hits of each position for 'G'],
                                                 [Hits of each position for 'A']]]
    And prints a dataframe with columns: sample, mir, offset, offsetNorm, counts, countsNorm, polarity
    in the outfile
    """
    dataframe = []
    dataframe.append("sample\tmir\toffset\toffsetNorm\tcounts\tcountsNorm\tpolarity")
    for ref in mirs.keys():
        """ For each reference name in mirs write the coverage of each of its positions """
        maximum = max(mirs[ref][2][0])
        for i in range(4):
            if max(mirs[ref][2][i]):
                maximum = max(mirs[ref][2][i])
        for pos in range(len(mirs[ref][2][1])):
            """ Compute coverage of each position and append to the dataframe a new line"""
            coverage = mirs[ref][2][0][pos]+mirs[ref][2][1][pos]+mirs[ref][2][2][pos]+mirs[ref][2][3][pos]
            dataframe.append("\t".join([sample, ref,
                                        str(pos+1), # offset
                                        str(float((pos+1))/len(mirs[ref][2][1])), # offsetNorm
                                        str(coverage), # count
                                        str(float(coverage)/float(maximum)) if maximum >0 else '0', # countNorm
                                        mirs[ref][1]])) # polarity
    try:
        out = open(outfile, 'w')
        out.write("\n".join(dataframe))
        out.write("\n")
        out.close()
    except IOError as e:
        sys.stderr.write("Error while writing file %s\n" % outfile)
        sys.stderr.write("I/O error(%s): %s\n" % (e.errno, e.strerror))

def write_counts(counts, outfile):
    """
    Takes a dictionary of counts[mir]=[counts] and prints it as a table of Gene    Counts
    """
    table = []
    table.append("Gene\tCounts")
    for gene in counts:
        table.append("\t".join([gene, str(counts[gene][0])]))
    try:
        out = open(outfile, 'w')
        out.write("\n".join(table))
        out.write("\n")
        out.close()
    except IOError as e:
        sys.stderr.write("Error while writing file %s\n" % outfile)
        sys.stderr.write("I/O error(%s): %s\n" % (e.errno, e.strerror))

def __main__():
    options = option_parsing()
    bamfile = pysam.AlignmentFile(options.alignment_file, 'rb', check_sq=False)
    pre_mirs = dict()
    mirs = dict()
    if options.pre_mirs:
        pre_mirs = get_count(bamfile)
        write_counts(pre_mirs, options.output_pre_mirs)
    if options.mirs:
        mirs = get_coverage(bamfile, options.gff_file)
        write_counts(mirs, options.output_mature_mirs)
    if options.lattice:
        write_dataframe(mirs, options.lattice, options.alignment_file)

if __main__():
    __main__()

