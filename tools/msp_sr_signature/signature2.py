import argparse
from collections import defaultdict

import numpy

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
        '--minscope', type=int,
        help="Minimum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--maxscope', type=int,
        help="Maximum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--output_h', action="store", type=str,
        help="h-signature dataframe")
    the_parser.add_argument(
        '--output_z', action="store", type=str,
        help="z-signature dataframe")
    the_parser.add_argument(
        '--graph', action='store', choices=["global", "by-item"],
        help="small RNA signature is computed either globally or\
        by item (global-lattice)")
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
        query_range = range (minquery, maxquery+1)
        target_range = range (mintarget, maxtarget+1)
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
                    Target_table[key[0]][coordinate] = Target_table[key[0]].get(
                        coordinate, 0) + 1
        return Query_table, Target_table

    def compute_signature_z(self, minquery, maxquery, mintarget, maxtarget,
                            scope, genome_wide=False, zscore="no"):
        Query_table, Target_table = signature_tables(minquery, maxquery,
                                                     mintarget, maxtarget)
        frequency_table = defaultdict(dict)
        for chrom in self.chromosomes:
            for overlap in scope:
                frequency_table[chrom][overlap] = 0
        for chrom in Query_table:
            for coord in Query_table[chrom]:
                for overlap in scope:
                    frequency_table[chrom][overlap] += min(
                        Query_table[chrom][coord],
                        Target_table[chrom].get(-coord -i +1, 0))
        # since we want the number of pairs, not the number or paired reads
        for chrom in frequency_table:
            for overlap in frequency_table[chrom]:
                frequency_table[chrom][overlap] /= 2
        # collapse if genome_wide = True
        for overlap in scope:
            accumulator = []
            for chrom in frequency_table:
                accumulator.append(frequency_table[chrom][overlap])
            frequency_table['all_chromosomes'][overlap] = sum(accumulator)
        # compute z-scores by chromosomes and for all-chromosomes
        for chrom in frequency_table:
            accumulator = []
            for overlap in frequency_table[chrom]:
                accumulator.append(frequency_table[chrom][overlap])
            z_mean = mean(accumulator)
            z_std = std(accumulator)
            if z_std == 0:
                for overlap in frequency_table[chrom]:
                    frequency_table[chrom][overlap] =
                        [frequency_table[chrom][overlap], 0]
            else:
                for overlap in frequency_table[chrom]:
                    frequency_table[chrom][overlap] =
                        [frequency_table[chrom][overlap],
                        (frequency_table[chrom][overlap]- z_mean)/z_std]
        return frequency_table

    def compute_signature_h_byitem(self, minquery, maxquery, mintarget,
                                   maxtarget, scope, genome_wide=False):
        Query_table, Target_table = signature_tables(minquery, maxquery,
                                                     mintarget, maxtarget)
        frequency_table = defaultdict(dict)
        for chrom in self.chromosomes:
            for overlap in scope:
                frequency_table[chrom][overlap] = 0
        for chrom in Query_table:
            Total_Query_Numb = 0
            for coord in Query_table[chrom]:
                Total_Query_Numb += Query_table[chrom][coord]
            for coord in Query_table[chrom]:
                local_table = dict ([(overlap, 0) for overlap in scope])
                number_of_targets = 0
                for overlap in scope:
                    local_table[overlap] += Query_table[chrom][coord] *
                        Target_table[chrom].get(-offset -i +1, 0)
                    number_of_targets +=
                        Target_table[chrom].get(-offset -i +1, 0)
                for overlap in scope:
                    try:
                        frequency_table[overlap] += (1. / number_of_targets
                            / Total_Query_Numb) * frequency_table[i]
                    except ZeroDivisionError :
                        continue
        return general_frequency_table      

    def format_table(self, table):
    

def main(input):
    F = open(output, 'w')
    I = open(input, 'read')
    mapobj = Map(input)
    something = mapobj.compute_signature_z(minquery, maxquery, mintarget,
                    maxtarget, scope, genome_wide=False, zscore="no")


if __name__ == "__main__":
    args = Parser()
    main(args., args., args., args.)