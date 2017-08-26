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

    def compute_signature_z(self, minquery, maxquery, mintarget, maxtarget,
                            scope, zscore="no"):
        Query_table, Target_table = self.signature_tables(minquery, maxquery,
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
                        Target_table[chrom].get(-coord - overlap + 1, 0))
        # since we want the number of pairs, not the number or paired reads
        #### WORK AT THIS, DO NOT LEAVE IT AS IS ####
        for chrom in frequency_table:
            for overlap in frequency_table[chrom]:
                frequency_table[chrom][overlap] /= 2
        # compute overlaps for all chromosomes merged
        for overlap in scope:
            accumulator = []
            for chrom in frequency_table:
                if chrom != 'all_chromosomes':
                    accumulator.append(frequency_table[chrom][overlap])
            frequency_table['all_chromosomes'][overlap] = sum(accumulator)
        return self.stringify_table(frequency_table)

    def compute_signature_h(self, minquery, maxquery, mintarget,
                                   maxtarget, scope):
        Query_table, Target_table = self.signature_tables(minquery, maxquery,
                                                     mintarget, maxtarget)
        frequency_table = defaultdict(dict)
        for chrom in self.chromosomes:
            for overlap in scope:
                frequency_table[chrom][overlap] = 0
        for chrom in Query_table:
            Total_Query_Numb = 0
            for coord in Query_table[chrom]:
                Total_Query_Numb += Query_table[chrom][coord]
                local_table = dict([(overlap, 0) for overlap in scope])
                number_of_targets = 0
                for overlap in scope:
                    local_table[overlap] += Query_table[chrom][coord] * \
                        Target_table[chrom].get(-coord - overlap + 1, 0)
                    number_of_targets += Target_table[chrom].get(
                        -coord - overlap + 1, 0)
                for overlap in scope:
                    try:
                        frequency_table[chrom][overlap] += \
                            local_table[overlap] / number_of_targets \
                            / float(Total_Query_Numb)
                    except ZeroDivisionError:
                        continue
        # compute overlap probabilities for all chromosomes merged
        general_frequency_table = dict([(overlap, 0) for overlap in scope])
        total_aligned_reads = 0
        for chrom in frequency_table:
            for overlap in frequency_table[chrom]:
                total_aligned_reads += self.bam_object.count(chrom)
        for chrom in frequency_table:
            for overlap in frequency_table[chrom]:
                try:
                    general_frequency_table[overlap] += \
                        frequency_table[chrom][overlap] / total_aligned_reads \
                        * self.bam_object.count(chrom)
                except ZeroDivisionError:
                    continue
        for overlap in general_frequency_table:
            frequency_table['all_chromosomes'][overlap] = \
                general_frequency_table[overlap]
        return self.stringify_table(frequency_table)

    def stringify_table(self, frequency_table):
        '''
        method both to compute z-score and to return a writable string
        '''
        tablestring = []
        for chrom in sorted(frequency_table):
            accumulator = []
            for overlap in frequency_table[chrom]:
                accumulator.append(frequency_table[chrom][overlap])
            z_mean = numpy.mean(accumulator)
            z_std = numpy.std(accumulator)
            if z_std == 0:
                for overlap in sorted(frequency_table[chrom]):
                    tablestring.append('%s\t%s\t%s\t%s\n' % (
                        chrom, str(overlap),
                        str(frequency_table[chrom][overlap]), str(0)))
            else:
                for overlap in sorted(frequency_table[chrom]):
                    tablestring.append('%s\t%s\t%s\t%s\n' % (
                        chrom, str(overlap),
                        str(frequency_table[chrom][overlap]),
                        str((frequency_table[chrom][overlap] - z_mean)/z_std)))
        return ''.join(tablestring)



def main(input, minquery, maxquery, mintarget, maxtarget, minscope, maxscope,
         output_h, output_z, genome_wide=False, zscore="no"):
    H = open(output_h, 'w')
    Z = open(output_z, 'w')
    mapobj = Map(input)
    scope = range(minscope, maxscope + 1)
    Z.write(mapobj.compute_signature_z(minquery, maxquery, mintarget,
            maxtarget, scope, zscore="no"))
    H.write(mapobj.compute_signature_h(minquery, maxquery, mintarget,
            maxtarget, scope))
    H.close()
    Z.close()


if __name__ == "__main__":
    args = Parser()
    main(args.input, args.minquery, args.maxquery, args.mintarget,
         args.maxtarget, args.minscope, args.maxscope, args.output_h,
         args.output_z)
