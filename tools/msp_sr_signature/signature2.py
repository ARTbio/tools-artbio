import argparse
from collections import defaultdict

import numpy

import pysam

def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="input alignment file")
    the_parser.add_argument(
        '--minquery', type=int, help="Minimum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxquery', type=int, help="Maximum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--mintarget', type=int, help="Minimum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxtarget', type=int, help="Maximum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--minscope', type=int, help="Minimum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--maxscope', type=int, help="Maximum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--outputOverlapDataframe', action="store", type=str, help="Overlap dataframe")
    the_parser.add_argument('--referenceGenome', action='store',
                            help="path to the bowtie-indexed or fasta reference")
    the_parser.add_argument('--extract_index', action='store_true',
                            help="specify if the reference is an indexed Bowtie reference")
    the_parser.add_argument('--graph', action='store', choices=[
                            "global", "lattice"], help="small RNA signature is computed either globally or by item (global-lattice)")
    the_parser.add_argument(
        '--rcode', type=str, help="R code to be passed to the python script")
    args = the_parser.parse_args()
    return args

class Map:

    def __init__(self, bam_file, sample):
        self.sample_name = sample
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
        Query_table = {}
        Target_table = {}
        frequency_table = {}
        for key in self.map_dict:
            for size in self.map_dict[key]:
                if size in query_range:
                    if key[2] = 'F':
                        coordinate = key[1]
                    else:
                        coordinate = -1 * key[1]
                    Query_table[key[0]][coordinate] = Query_table[key[0]].get(
                        coordinate, 0) + 1
                if size in target_range:
                    if key[2] = 'F':
                        coordinate = key[1]
                    else:
                        coordinate = -1 * key[1]
                    Target_table[key[0]][coordinate] = Target_table[key[0]].get(
                        coordinate, 0) + 1
        return Query_table, Target_table

    def compute_signature_z(self, minquery, maxquery, mintarget, maxtarget,
                            scope, genome_wide=False, zscore="no"):
        for chrom in self.chromosomes:
            for overlap in scope:
                frequency_table[chrom][overlap] = 0
            for key in Query_table:
                for i in scope:
                    if key[2]== 'F':
                        frequency_table[chrom][i] += min(Query_table[key],
                            Target_table.get((chrom, key[1] + i - 1, 'R'),
                                              0))
                    else:
                        frequency_table[chrom][i] += min(Query_table[key],
                            Target_table.get((chrom, key[1] - i + 1, 'F'),
                                              0))
        # since we want the number of pairs, not the number or paired reads
        for chrom in frequency_table:
            for i in frequency_table[chrom]:
                frequency_table[chrom][i] /= 2
        # collapse if genome_wide = True
        for i in scope:
            accumulator = []
            for chrom in frequency_table:
                accumulator.append(frequency_table[chrom][i])
            frequency_table['all_chromosomes'][i] = sum(accumulator)
        # compute z-scores by chromosomes and for all-chromosomes
        for chrom in frequency_table:
            accumulator = []
            for i in frequency_table[chrom]:
                accumulator.append(frequency_table[chrom][i])
            z_mean = mean(accumulator)
            z_std = std(accumulator)
            if z_std == 0:
                for i in frequency_table[chrom]:
                    frequency_table[chrom][i] = [frequency_table[chrom][i], 0]
            else:
                for i in frequency_table[chrom]:
                    frequency_table[chrom][i] = [frequency_table[chrom][i],
                        (frequency_table[chrom][i]- z_mean)/z_std]
        return frequency_table

    def compute_signature_h(self, minquery, maxquery, mintarget, maxtarget,
                            scope, genome_wide=False):
        query_range = range (minquery, maxquery+1)
        target_range = range (mintarget, maxtarget+1)
        Query_table = {}
        Target_table = {}
        Total_Query_Numb = 0
        frequency_table = {}
        for chrom in self.chromosomes:
            for overlap in scope:
                frequency_table[chrom][overlap] = 0
        for chrom in self.chromosomes:
            for key in self.readDict:
                for size in self.readDict[key]:
                    if size in query_range:
                        Query_table[key] = Query_table.get(key, 0) + 1
                    if size in target_range:
                        Target_table[key] = Target_table.get(key, 0) + 1
    for offset in Query_table:
      frequency_table = dict ([(i,0) for i in scope])
      number_of_targets = 0
      for i in scope:
        frequency_table[i] += Query_table[offset] *  Target_table.get(-offset -i +1, 0)
        number_of_targets += Target_table.get(-offset -i +1, 0)
      for i in scope:
        try:
          general_frequency_table[i] += (1. / number_of_targets / Total_Query_Numb) * frequency_table[i]
        except ZeroDivisionError :
          continue
    return general_frequency_table      


def main():
    F = open(output, 'w')
    I = open(input, 'read')
    mapobj = Map(input, sample)
    something = mapobj.compute_signature_z(minquery, maxquery, mintarget,
                    maxtarget, scope, genome_wide=False, zscore="no")


if __name__ == "__main__":
    args = Parser()
    # if identical sample names
    if len(set(args.sample_names)) != len(args.sample_names):
        args.sample_names = [name + '_' + str(i) for
                             i, name in enumerate(args.sample_names)]
    main(args., args., args., args.)