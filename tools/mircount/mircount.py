#!/usr/bin/env python
# python parser module for pre-mir and mature miRNAs, guided by mirbase.org GFF3
# version 0.0.9 (1-6-2014)
# Usage MirParser.py  <1:index source> <2:extraction directive> <3:output pre-mir> <4: output mature miRs> <5:mirbase GFF3>
#                     <6:pathToLatticeDataframe or "dummy_dataframe_path"> <7:Rcode or "dummy_plotCode"> <8:latticePDF or "dummy_latticePDF">
#                     <9:10:11 filePath:FileExt:FileLabel> <.. ad  lib>

import sys
import subprocess
import sys
import logging
import optparse
import time

#from smRtools import HandleSmRNAwindows
""" Log parsing """
LOG_FORMAT = '%(asctime)s|%(levelname)-8s|%(message)s'
LOG_DATEFMT = '%Y-%m-%d %H:%M:%S'
LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


def getHeaders(fasta_file, logger):
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
        if polarity == "R":
            self.aligned_reads[offset].append(size)
        else:
            self.aligned_reads[-(offset + size - 1)].append(size)
        self.aligned_reads_count += 1

    def count_reads(self, upstream_coord, downstream_coord, polarity):
        """
        Takes upstream and downstream coordinates as well as the polarity of
        the pre-mir and counts the number of reads overlaping it.
        """
        overlaping = 0
        if polarity == "F":
            """
            If the polarity is + then check for each position
            between upstream and downstream coord if a read has a 5' matching
            """
            for offset in xrange(upstream_coord, downstream_coord + 1)
                if self.aligned_reads.has_key(offset):
                    for read in self.aligned_reads[offset]:
                        overlaping += 1
        else:
            """
            If the polarity is + then check for each position
            between upstream and downstream coord if a read has a 3' matching
            """
            for offset in xrange(upstream_coord, downstream_coord + 1)
                if self.aligned_reads.has_key(-offset):
                    for read in self.aligned_reads[-offset]:
                        overlaping += 1
        return overlaping

def __main__():
    """ Get options """
    parser = optparse.OptionParser(description='Count miRNAs')
    parser.add_option('--ref_genome', type='string', dest='ref_genome_file', help='Fastq file of clipped sequence reads', metavar='FILE')
    parser.add_option('--alignment', type='string', dest='alignment_file', help='Alignment tabular or sam file')
    parser.add_option('--alignment_format', choices=['tabular', 'sam', 'bam'], dest='alignment_format', help='Alignement format (tabular, sam or bam). [default=tabular]', default='tabular')
    parser.add_option('--pre_mirs_output', type='string', dest='output_pre_mirs', help='GFF3 describing the mature miRs', metavar='FILE')
    parser.add_option('--mature_mirs_output', type='string', dest='output_mature_mirs', help='Reference genome', metavar='FILE')
    parser.add_option('--gff', type='string', dest='gff_file', help=' GFF3 describing both pre-mirs and mature mirs (to be discussed)', metavar='FILE')
    parser.add_option('--lattice', type='string', dest='lattice')
    parser.add_option('--loglevel', choices=LOG_LEVELS, default='INFO', help='logging level (default: INFO)')
    parser.add_option('-l', '--logfile', help='log file (default=stderr)')
    (options, args) = parser.parse_args()
    
    """ Check if the options were correctly passed """
    if len(args) > 0:
        parser.error('Wrong number of arguments')
    if (not options.ref_genome_file or not options.alignment_file or 
        not options.gff_file:# or not options._file or not options.ref_genome_gff_file):
        parser.error('Missing file')
    
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
    headers = list()
    hit_store = dict()
    polarity = ''
    gene = ''
    offset = None
    size = None
    """ Initialize HitContainerStore """
    for header in headers:
        hit_store[header] = HitContainer(header)
    """!!! Get reference headers and read alignment file !!!"""
    try:
        """ Get fasta headers """
        headers = getHeaders(ref_genome_file, logger)
        """ If the pre-mirs are needed """
        if options.output_pre_mirs :
            if options.alignment_format != 'tabular':
                import pysam
                samfile = None
                """ Read file with pysam """
                if options.alignment_format == "bam" :
                    samfile = pysam.AlignmentFile(options.alignment_file, 'rb')
                else:
                    samfile = pysam.AlignmentFile(options.alignment_file, 'r')
                for read in samfile:
                    """
                    Get:
                     - read polarity
                     - reference name
                     - read offset (corrected 0-based coordinates to 1-based)
                     - read length
                    """
                    polarity = "F" if not read.is_reverse else polarity = "R"
                    gene = read.reference_name
                    offset = read.pos + 1
                    size = read.qlen
                    hit_store[gene].add_read(offset, polarity, size)
            else:
                fh = open(options.alignment_file, 'r')
                for line in fh.readlines():
                    fields = line.split()
                    polarity = fields[1]
                    gene = fields[2]
                    offset = int(fields[3])
                    size = len(fields[4])
                    hit_store[gene].add_read(offset, polarity, size)
                fh.close()
    except IOError as e:
        raise e
    """ Read GFF file and count hits """
    

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
