import argparse
from collections import defaultdict

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="bam alignment file")
    the_parser.add_argument(
        '--minquery', type=int,
        help="Minimum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxquery', type=int,
        help="Maximum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--mintarget', type=int,
        help="Minimum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxtarget', type=int,
        help="Maximum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--overlap', type=int,
        help="Overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--output', action="store", type=str,
        help="Pairable sequences")
    args = the_parser.parse_args()
    return args


class Map:

    def __init__(self, bam_file):
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object)

    def create_map(self, bam_object):
        '''
        Returns a map_dictionary {(chromosome,read_position,polarity):
                                                    [read_length, ...]}
        '''
        map_dictionary = defaultdict(list)
        # get empty value for start and end of each chromosome
        for chrom in self.chromosomes:
            map_dictionary[(chrom, 1, 'F')] = []
            map_dictionary[(chrom, self.chromosomes[chrom], 'F')] = []
        for chrom in self.chromosomes:
            for read in bam_object.fetch(chrom):
                positions = read.positions  # a list of covered positions
                if read.is_reverse:
                    map_dictionary[(chrom, positions[-1]+1,
                                    'R')].append(read.query_alignment_length)
                else:
                    map_dictionary[(chrom, positions[0]+1,
                                    'F')].append(read.query_alignment_length)
        return map_dictionary

    def signature_tables(self, minquery, maxquery, mintarget, maxtarget):
        query_range = range(minquery, maxquery + 1)
        target_range = range(mintarget, maxtarget + 1)
        Query_table = defaultdict(dict)
        Target_table = defaultdict(dict)
        for key in self.map_dict:
            for size in self.map_dict[key]:
                if size in query_range or size in target_range:
                    if key[2] == 'F':
                        coordinate = key[1]
                    else:
                        coordinate = -key[1]
                if size in query_range:
                    Query_table[key[0]][coordinate] = Query_table[key[0]].get(
                        coordinate, 0) + 1
                if size in target_range:
                    Target_table[key[0]][coordinate] = \
                        Target_table[key[0]].get(coordinate, 0) + 1
        return Query_table, Target_table

    def search_overlaps(self, minquery, maxquery, mintarget, maxtarget,
                        overlap=10):
        Query_table, Target_table = self.signature_tables(minquery, maxquery,
                                                          mintarget, maxtarget)
        overlap_groups = defaultdict(list)
        for chrom in Query_table:
            for coord in Query_table[chrom]:
                if Target_table[chrom].get(-coord - overlap + 1, 0):
                    overlap_groups[chrom].append(coord)
        return overlap_groups

    def feed_overlaps(self, overlap_groups, minquery, output, overlap=10):
        F = open(output, 'w')
        for chrom in sorted(overlap_groups):
            for pos in sorted(overlap_groups[chrom]):
                if pos > 0: # read are forward
                    reads = self.bam_object.fetch(chrom, start=pos-1,
                                                  end=pos-1+overlap-1)
                    for read in reads:
                        positions = read.positions
                        if pos-1 == positions[0] and \
                            read.query_alignment_length >= minquery:
                                F.write('%s\t%s\t%s\n' % (
                                    chrom, pos, read.query_sequence))
                else: # reads are reverse
                    reads = self.bam_object.fetch(chrom,
                                                  start=-pos-1-overlap+1,
                                                  end=-pos-1)
                    for read in reads:
                        positions = read.positions
                        if -pos-1 == positions[-1] and \
                            read.query_alignment_length >= minquery:
                                F.write('%s\t%s\t%s\n' % (chrom, pos,
                                                      read.query_sequence))
        return

def main(input, minquery, maxquery, mintarget, maxtarget, output, overlap=10):
    mapobj = Map(input)
    mapobj.feed_overlaps(mapobj.search_overlaps(minquery, maxquery,
                                                mintarget, maxtarget,
                                                overlap), minquery, output)

if __name__ == "__main__":
    args = Parser()
    main(args.input, args.minquery, args.maxquery, args.mintarget,
         args.maxtarget, args.output)
