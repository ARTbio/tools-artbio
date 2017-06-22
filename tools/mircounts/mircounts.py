#!/usr/bin/env python
# python parser module for pre-mir and mature miRNAs, guided by mirbase.org GFF3
# version 0.0.9 (1-6-2014)
# Usage MirParser.py  <1:index source> <2:extraction directive> <3:output pre-mir> <4: output mature miRs> <5:mirbase GFF3>
#                     <6:pathToLatticeDataframe or "dummy_dataframe_path"> <7:Rcode or "dummy_plotCode"> <8:latticePDF or "dummy_latticePDF">
#                     <9:10:11 filePath:FileExt:FileLabel> <.. ad  lib>

import sys
import logging
import optparse
import time
from collections import defaultdict
#from smRtools import HandleSmRNAwindows
""" Log parsing """
LOG_FORMAT = '%(asctime)s|%(levelname)-8s|%(message)s'
LOG_DATEFMT = '%Y-%m-%d %H:%M:%S'
LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

def get_headers_from_fasta(fasta_file, logger):
    """
    Takes a fasta file and returns a list of fasta headers.
    """
    headers = list()
    try:
        fh = open(fasta_file, 'r')
        for line in fh.readlines():
            if line.startswith('>', 0, 1):
                """ Bowtie splits headers """
                headers.append(line[1:-1].split()[0])
        fh.close()
    except IOError as e:
        logger.error("I/O error(%s): %s" % (e.errno, e.strerror))
        raise e
    return headers

def read_bam_sam(alignment_file, alignment_format):
    """
    Takes an alignment file and its format (must be either bam or sam)
    and returns a HitContainer dictionary
    """
    import pysam
    hit_store = dict()
    gene = None
    offset = None
    size = None
    samfile = None
    samfile = pysam.AlignmentFile(alignment_file, 'r')
    """ Initialize dictionary """
    for reference in samfile.header['SQ']:
        hit_store[reference['SN']] = HitContainer(reference['SN'])
    for read in samfile:
        """
        Get:
         - read polarity
         - reference name
         - read offset (corrected 0-based coordinates to 1-based)
         - read length
        """
        polarity = "forward" if (not read.is_reverse) else "reverse"
        if not read.is_unmapped:
            gene = read.reference_name
            offset = read.pos + 1
            size = read.qlen
            """ If the gene hasn't been read yet add it to store """
            hit_store[gene].add_read(offset, polarity, size)
    return hit_store

def read_tabular(alignment_file):
    """
    Reads tabular alignment file and returns HitContainer dictionary
    """
    hit_store = dict()
    fields = None
    polarity = None
    gene = None
    offset = None
    size = None
    try:
        fh = open(alignment_file, 'r')
        for line in fh.readlines():
            fields = line.split()
            polarity = fields[1]
            gene = fields[2]
            offset = int(fields[3])
            size = len(fields[4])
            try:
                hit_store[gene].add_read(offset, polarity, size)
            except:
                hit_store[gene] = HitContainer(gene)
                hit_store[gene].add_read(offset, polarity, size)
        fh.close()
    except IOError as e:
        raise e
    return hit_store

class HitContainer:
    """
    This is a class made to store each read that hit a particular fasta sequence during alignment
    """
    def __init__(self, header):
        self.fasta_seq = header
        self.aligned_reads_count = 0
        self.aligned_reads = defaultdict(list)

    def add_read(self, offset, polarity, size):
        """
        Add read to the self.aligned_reads dictionary as self.aligned_reads[offset] = [size]
        """
        if polarity == "forward":
            self.aligned_reads[offset].append(size)
        else:
            self.aligned_reads[-(offset + size - 1)].append(size)
        self.aligned_reads_count += 1

    def count_reads(self, upstream_coord, downstream_coord, polarity):
        """
        Takes upstream and downstream coordinates as well as the polarity of
        the mir and counts the number of reads overlaping it.
        """
        overlaping = 0
        if polarity == "forward":
            """
            If the polarity is + then check for each position
            between upstream and downstream coord if a read has a 5' matching
            """
            for offset in xrange(upstream_coord, downstream_coord + 1):
                if self.aligned_reads.has_key(offset):
                    for read in self.aligned_reads[offset]:
                        overlaping += 1
        else:
            """
            If the polarity is + then check for each position
            between upstream and downstream coord if a read has a 3' matching
            """
            for offset in xrange(upstream_coord, downstream_coord + 1):
                if self.aligned_reads.has_key(-offset):
                    for read in self.aligned_reads[-offset]:
                        overlaping += 1
        return overlaping

def __main__():
    """ Get options """
    parser = optparse.OptionParser(description='Count miRNAs')
    parser.add_option('--loglevel', choices=LOG_LEVELS, default='INFO', help='logging level (default: INFO)')
    parser.add_option('-p', '--pre_mirs', action='store_true', dest='pre_mirs', help='Count pre-miRNAs', default=False)
    parser.add_option('-m', '--mirs', action='store_true', dest='mirs', help='Count mature miRNAs', default=False)
    inputs = optparse.OptionGroup(parser, 'Inputs')
    inputs.add_option('--alignment', type='string', dest='alignment_file', help='Alignment tabular or sam file')
    inputs.add_option('--alignment_format', choices=['tabular', 'sam', 'bam'], dest='alignment_format', help='Alignement format (tabular, sam or bam). [default=tabular]', default='tabular')
    inputs.add_option('--gff', type='string', dest='gff_file', help=' GFF3 describing both pre-mirs and mature mirs (to be discussed)', metavar='FILE')
    parser.add_option_group(inputs)
    outputs = optparse.OptionGroup(parser, 'Outputs')
    outputs.add_option('--pre_mirs_output', type='string', dest='output_pre_mirs', help='GFF3 describing the mature miRs', metavar='FILE', default='output_pre_mirs_count.tab')
    outputs.add_option('--mature_mirs_output', type='string', dest='output_mature_mirs', help='Reference genome', metavar='FILE', default='output_mirs_count.tab')
    outputs.add_option('--lattice', type='string', dest='lattice')
    outputs.add_option('-l', '--logfile', help='log file (default=stderr)')
    parser.add_option_group(outputs)
    (options, args) = parser.parse_args()
    """ Check if the options were correctly passed """
    if len(args) > 0:
        parser.error('Wrong number of arguments')
    if (not options.alignment_file or not options.gff_file):
        parser.error('Missing file')
    if (options.lattice and not options.mirs):
        parse.error("Can't output lattice dataframe if the option '-m' is not set")
    """ Set up the logger """
    log_level = getattr(logging, options.loglevel)
    kwargs = {'format': LOG_FORMAT,
              'datefmt': LOG_DATEFMT,
              'level': log_level}
    if options.logfile:
        kwargs['filename'] = options.logfile
    logging.basicConfig(**kwargs)
    logger = logging.getLogger('MirCount')
    
    """ Declare Variables """
    hit_store = dict()
    """ Read alignment file """
    if options.alignment_format != 'tabular':
        hit_store = read_bam_sam(options.alignment_file, options.alignment_format)
    else:
        hit_store = read_tabular(options.alignment_file)
        logger.info("Number of keys %s" % len(hit_store.keys()))
    """ If the pre-mirs are asked """
    if options.pre_mirs:
        text = list()
        out_pre_mirs = options.output_pre_mirs
        try:
            """ header """
            text.append("Gene\tCount")
            for gene in hit_store.keys():
                text.append("\t".join([gene,str(hit_store[gene].aligned_reads_count)]))
            fh = open(out_pre_mirs, 'w')
            fh.write("\n".join(text))
            fh.close()
            logger.info("Finished writing %s" % out_pre_mirs)
        except IOError as e:
            logger.error("I/O error(%s): %s" % (e.errno, e.strerror))
    """ Read GFF file and count hits """
    if options.mirs:
        out_mirs = options.output_mature_mirs
        text = list()
        if options.lattice:
            lattice_dataframe = list()
            lattice_dataframe.append("sample\tmir\toffset\toffsetNorm\tcounts\tcountsNorm\tpolarity")
            coverage = None
        try:
            """ Open GFF and Output file """
            gff = open(options.gff_file, 'r')
            text.append("Gene\tCount")
            for line in gff.readlines():
                """ For each line if it doesn't start with '#' split it and get fields """
                if line[0] != '#':
                    gff_fields = line[:-1].split("\t")
                    if gff_fields[2] == 'miRNA':
                        chrom = gff_fields[0]
                        if hit_store.has_key(chrom):
                            item_upstream_coord = int(gff_fields[3])
                            item_downstream_coord = int(gff_fields[4])
                            if gff_fields[6] == '+':
                                item_polarity = 'forward'
                                if options.lattice:
                                    """
                                    If the lattice dataframe is needed as output
                                    we need to compute the coverage off each GFF item
                                     - initialize coverage dictionary with each position of the item as key and 0 as value
                                     - for each offset hit by a read for this item:
                                       - get the size of the reads that hit said offset
                                       - set +1 for each position covered by each read in coverage dictionary
                                    """
                                    coverage = dict([(i,0) for i in xrange(1, item_downstream_coord - item_upstream_coord+1)])
                                    for offset in hit_store[chrom].aligned_reads.keys():
                                        if (offset > 0) and (coverage.has_key(offset - item_upstream_coord+1)):
                                            for read in hit_store[chrom].aligned_reads[offset]:
                                                it = 0
                                                while it < read:
                                                    if coverage.has_key(offset - item_upstream_coord+1+it):
                                                        coverage[offset - item_upstream_coord+1+it] += 1
                                                        it += 1
                                                    else:
                                                        it = read
                            else:
                                item_polarity = 'reverse'
                                if options.lattice:
                                    """
                                    If the lattice dataframe is needed as output
                                    we need to compute the coverage off each GFF item
                                     - initialize coverage dictionary with each position of the item as key and 0 as value
                                     - for each offset hit by a read for this item:
                                       - get the size of the reads that hit said offset
                                       - set +1 for each position covered by each read in coverage dictionary
                                    """
                                    coverage = dict([(i,0) for i in xrange(1, item_downstream_coord-item_upstream_coord+1)])
                                    for offset in hit_store[chrom].aligned_reads.keys():
                                        if (offset < 0) and (coverage.has_key(-offset - item_upstream_coord+1)):
                                            for read in hit_store[chrom].aligned_reads[offset]:
                                                it = 0
                                                while it < read:
                                                    if coverage.has_key(-offset - item_upstream_coord - it):
                                                        coverage[-offset - item_upstream_coord - it] += 1
                                                        it += 1
                                                    else:
                                                        it = read
                            """ count reads and write table """
                            count = hit_store[chrom].count_reads(item_upstream_coord,
                                                                 item_downstream_coord, item_polarity)
                            text.append("\t".join([chrom, str(count)]))
                            if options.lattice:
                                """ Create a line in the lattice_dataframe for each position of the item """
                                maximum = max(coverage.values())
                                lattice_lines = list()
                                for pos in coverage.keys():
                                    lattice_lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (options.alignment_file,
                                                                                         chrom, pos,
                                                                                         float(pos)/(item_downstream_coord - item_upstream_coord +1),
                                                                                         coverage[pos],
                                                                                         (float(coverage[pos])/maximum) if maximum > 0 else 0,
                                                                                         item_polarity))
                                lattice_dataframe.append("\n".join(lattice_lines))
            gff.close()
            fh_out_mirs = open(out_mirs, 'w')
            fh_out_mirs.write("\n".join(sorted(text)))
            fh_out_mirs.close()
            logger.info("Finished writing %s" % out_mirs)
            if options.lattice:
                """ Print lattice dataframe if asked for """
                fh_out_lattice = open(options.lattice, 'w')
                fh_out_lattice.write("\n".join(lattice_dataframe))
                fh_out_lattice.close()
                logger.info("Finished writing %s" % options.lattice)
        except IOError as e:
            logger.error("I/O error(%s): %s" % (e.errno, e.strerror))
        except KeyError as e:
            logger.error("The first column of the GFF and the headers of the reference genome of the alignment are not the same")
            logger.error("We caught %s in the GFF3 and an example of reference header would be %s" % (e, hit_store.keys()[0]))

def __test__():
    aln = sys.argv[1]
    aln_fmt = sys.argv[2]
    print aln,aln_fmt
    hit_store = read_bam_sam(aln,aln_fmt)
    print "NB genes:\t",len(hit_store.keys())
    print "Gene\tCount"
    for gene in hit_store.keys():
        print gene+"\t"+str(hit_store[gene].aligned_reads_count)

if __main__():
    __main__()
"""
IndexSource = sys.argv[1]
ExtractionDirective = sys.argv[2]
if ExtractionDirective == "--do_not_extract_index":
    genomeRefFormat = "fastaSource"
elif ExtractionDirective == "--extract_index":
    genomeRefFormat = "bowtieIndex"
OutputPre_mirs = sys.argv[3]
OutputMature_Mirs = sys.argv[4]
GFF3_file = sys.argv[5]
lattice = sys.argv[6]
Rcode = sys.argv[7]
latticePDF = sys.argv[8]
Triplets = [sys.argv[9:][i:i + 3] for i in xrange(0, len(sys.argv[9:]), 3)]
MasterListOfGenomes = {}

for [filePath, FileExt, FileLabel] in Triplets:
    print FileLabel
    MasterListOfGenomes[FileLabel] = HandleSmRNAwindows(alignmentFile=filePath,
                                                        alignmentFileFormat=FileExt,
                                                        genomeRefFile=IndexSource,
                                                        genomeRefFormat=genomeRefFormat,
                                                        biosample=FileLabel)

header = ["gene"]
for [filePath, FileExt, FileLabel] in Triplets:
    header.append(FileLabel)

hit_table = ["\t".join(header)]  # table header: gene, sample1, sample2, sample3, etc. separated by tabulation

# read GFF3 to subinstantiate
gff3 = open(GFF3_file, "r")
lattice_dataframe = []
for line in gff3:
    if line[0] == "#":
        continue
    gff_fields = line[:-1].split("\t")
    chrom = gff_fields[0]
    gff_name = gff_fields[-1].split("Name=")[-1].split(";")[0]  # to isolate the GFF Name
    item_upstream_coordinate = int(gff_fields[3])
    item_downstream_coordinate = int(gff_fields[4])
    if gff_fields[6] == "+":
        item_polarity = "forward"
    else:
        item_polarity = "reverse"
    item_line = [gff_name]
    for sample in header[1:]:
        count = MasterListOfGenomes[sample].instanceDict[chrom].readcount(upstream_coord=item_upstream_coordinate,
                                                                          downstream_coord=item_downstream_coordinate,
                                                                          polarity=item_polarity)
        item_line.append(str(count))
        # subtreatement for lattice
        if lattice != "dummy_dataframe_path":
            if ("5p" not in gff_name) and ("3p" not in gff_name):
                lattice_dataframe.append(MasterListOfGenomes[sample].instanceDict[chrom].readcoverage(
                    upstream_coord=item_upstream_coordinate,
                    downstream_coord=item_downstream_coordinate,
                    windowName=gff_name + "_" + sample))
        # end of subtreatement for lattice
    hit_table.append("\t".join(item_line))
gff3.close()

Fpremirs = open(OutputPre_mirs, "w")
print >> Fpremirs, hit_table[0]
finalPreList = [i for i in sorted(hit_table[1:]) if ("5p" not in i) and ("3p" not in i)]
print >> Fpremirs, "\n".join(finalPreList)
Fpremirs.close()

Fmaturemires = open(OutputMature_Mirs, "w")
print >> Fmaturemires, hit_table[0]
finalMatureList = [i for i in sorted(hit_table[1:]) if ("5p" in i) or ("3p" in i)]
print >> Fmaturemires, "\n".join(finalMatureList)
Fmaturemires.close()

if lattice != "dummy_dataframe_path":
    Flattice = open(lattice, "w")
    print >> Flattice, "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("sample",
                                                       "mir",
                                                       "offset",
                                                       "offsetNorm",
                                                       "counts",
                                                       "countsNorm",
                                                       "polarity")
    print >> Flattice, "\n".join(lattice_dataframe)
    Flattice.close()
    R_command = "Rscript " + Rcode
    process = subprocess.Popen(R_command.split())
    process.wait()
"""
