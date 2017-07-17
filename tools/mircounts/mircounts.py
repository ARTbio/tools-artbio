#!/usr/bin/env python
import optparse
from collections import defaultdict
import pysam
import sys

def option_parsing():
    parser = optparse.OptionParser(description='Count miRNAs')
    parser.add_option('-p', '--pre_mirs', action='store_true',
                      dest='pre_mirs', help='Count pre-miRNAs', default=False)
    parser.add_option('-m', '--mirs', action='store_true', dest='mirs',
                      help='Count mature miRNAs', default=False)
    inputs = optparse.OptionGroup(parser, 'Inputs')
    inputs.add_option('--alignment', type='string', dest='alignment_file',
                      help='Alignment tabular or sam file')
    inputs.add_option('--gff', type='string', dest='gff_file',
                      help='GFF3 describing both pre-miRNAs and mature miRNAs',
                      metavar='FILE')
    inputs.add_option('--quality_threshold', type='int', dest='quality_threshold',
                      help='Quality threshold for coverage (default=10)',
                      default=10)
    parser.add_option_group(inputs)
    outputs = optparse.OptionGroup(parser, 'Outputs')
    outputs.add_option('--pre_mirs_output', type='string', dest='output_pre_mirs',
                       help='pre-miRNA count file', metavar='FILE',
                       default='output_pre_mirs_count.tab')
    outputs.add_option('--mature_mirs_output', type='string',
                       dest='output_mature_mirs',
                       help='Mature miRNA counts file', metavar='FILE',
                       default='output_mirs_count.tab')
    outputs.add_option('--lattice', type='string', dest='lattice',
                       help='Output file for the lattice dataframe. Default: None')
    parser.add_option_group(outputs)
    (options, args) = parser.parse_args()
    """ Check if the options were correctly passed """
    if len(args) > 0:
        parser.error('Please use the options flags to pass your arguments')
    if (not options.alignment_file or not options.gff_file):
        parser.error("Missing file. Both '--alignment' and '--gff' files are needed")
    return options

def get_pre_mir_counts(bamfile):
    """
    Takes a AlignmentFile object and returns a dictionary of counts for reads
    aligning with pre_mirs (as keys)
    """
    count = dict()
    for ref_name in bamfile.references:
        count[ref_name] = bamfile.count(reference = ref_name)
    return count
    
def get_pre_mir_coverage(bamfile, quality=10)
    """
    Takes a AlignmentFile object and returns a dictionary of lists
    of coverage along the coordinates of pre_mirs (as keys)
    """
    coverage = dict()
    for ref_name, ref_len in zip(bamfile.references, bamfile.lengths):
        coverage[ref_name]= bamfile.count_coverage(reference = ref_name,
	                                               start=0, end = ref_len,
	                                               quality_threshold =quality)]
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
    counts = defaultdict(int)
    refs = bamfile.references
    for line in open(gff_file, 'r'):
        if line[0] != '#':
            gff_fields = line[:-1].split("\t")
            if gff_fields[2] == 'miRNA':
                mir_name = gff_fields[0]
                premir_name = gff_fields[8].split('=')[1]
                mir_start = int(gff_fields[3])
                mir_end = int(gff_fields[4])
                count = 0
                # GFF is 1-based, pysam is 0-based.
                counts[mir_name] = bamfile.count(reference = mir_name,
                                                 start = mir_start-1,
                                                 end = mir_end-1)
    return counts

def write_dataframe_coverage(countdict, outfile):
    """
    Takes a dict[pre_mir reference name] = [coverage list]
    and writes a dataframe with columns:
    <gene_type name>, offset, normoffset, counts and normcounts
    in the outfile
    """
    F = open(outfile, 'w')
    F.write('Pre-Mir\tOffset\tNorm_offset\tCount\tNorm_count')
    for ref in sorted(countdict):
        """ For each reference name in mirs write the coverage of each of its positions """
        max = max(countdict[ref])
        reference_length = len(countdict[ref])
        for pos, c in enumerate(countdict[ref]):
            """ Compute and write value for each reference position"""
            F.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ref,
                                str(pos + 1),
                                str(float(pos+1)/reference_length),
                                str(c),
                                str(float(c)/max) if max != 0 else '0'))

def write_counts(countdict, outfile, gene_type='mir'):
    """
    Takes a dict[<gene_type name>]=count and 
    writes a count table
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
        pre_mirs = get_pre_mir_counts(bamfile, options.quality_threshold)
        write_counts(pre_mirs, options.output_pre_mirs)
        if options.lattice:
            write_dataframe(pre_mirs, options.lattice, options.sample_name)
    if options.mirs:
        mirs = get_mir_counts(bamfile, options.gff_file, options.quality_threshold)
        write_counts(mirs, options.output_mature_mirs)

if __main__():
    __main__()

