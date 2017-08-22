import argparse
from collections import defaultdict

import numpy

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--input', dest='input', required=True,
                            nargs='+', help='input BAM files')
    the_parser.add_argument('--sample_name', dest='sample_name',
                            required=True, nargs='+', help='sample name')
    the_parser.add_argument('--output', action='store',
                            type=str, help='output tabular file')
    the_parser.add_argument('-S', '--sizes', action='store',
                            help='use to output read sizes dataframe')
    args = the_parser.parse_args()
    return args


class Map:

    def __init__(self, bam_file, sample, computeSize=False):
        self.sample_name = sample
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object)
        self.max = self.compute_max(self.map_dict)
        self.mean = self.compute_mean(self.map_dict)
        self.median = self.compute_median(self.map_dict)
        self.coverage = self.compute_coverage(self.map_dict)
        if computeSize:
            self.size = self.compute_size(self.map_dict)

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
                for pos in positions:
                    if not map_dictionary[(chrom, pos+1, 'F')]:
                        map_dictionary[(chrom, pos+1, 'F')] = []
                if read.is_reverse:
                    map_dictionary[(chrom, positions[-1]+1,
                                    'R')].append(read.query_alignment_length)
                else:
                    map_dictionary[(chrom, positions[0]+1,
                                    'F')].append(read.query_alignment_length)
        return map_dictionary

    def compute_max(self, map_dictionary):
        '''
        takes a map_dictionary as input and returns
        a max_dictionary {(chromosome,read_position,polarity):
                              max_of_number_of_read_at_any_position}
        '''
        merge_keylist = [(i[0], 0) for i in map_dictionary.keys()]
        max_dictionary = dict(merge_keylist)
        for key in map_dictionary:
            if len(map_dictionary[key]) > max_dictionary[key[0]]:
                max_dictionary[key[0]] = len(map_dictionary[key])
        return max_dictionary

    def compute_mean(self, map_dictionary):
        '''
        takes a map_dictionary as input and returns
        a mean_dictionary {(chromosome,read_position,polarity):
                                                mean_value_of_reads}
        '''
        mean_dictionary = dict()
        for key in map_dictionary:
            if len(map_dictionary[key]) == 0:
                mean_dictionary[key] = 0
            else:
                mean_dictionary[key] = round(numpy.mean(map_dictionary[key]),
                                             1)
        return mean_dictionary

    def compute_median(self, map_dictionary):
        '''
        takes a map_dictionary as input and returns
        a mean_dictionary {(chromosome,read_position,polarity):
                                                    mean_value_of_reads}
        '''
        median_dictionary = dict()
        for key in map_dictionary:
            if len(map_dictionary[key]) == 0:
                median_dictionary[key] = 0
            else:
                median_dictionary[key] = numpy.median(map_dictionary[key])
        return median_dictionary

    def compute_coverage(self, map_dictionary, quality=10):
        '''
        takes a map_dictionary as input and returns
        a coverage_dictionary {(chromosome,read_position,polarity):
                                                coverage}
        '''
        coverage_dictionary = dict()
        for chrom in self.chromosomes:
            coverage_dictionary[(chrom, 1, 'F')] = 0
            coverage_dictionary[(chrom, self.chromosomes[chrom], 'F')] = 0
        for key in map_dictionary:
            coverage = self.bam_object.count_coverage(
                                                reference=key[0],
                                                start=key[1]-1,
                                                end=key[1],
                                                quality_threshold=quality)
            """ Add the 4 coverage values """
            coverage = [sum(x) for x in zip(*coverage)]
            coverage_dictionary[key] = coverage[0]
            # coverage_dictionary[(key[0], key[1], 'R')] = coverage
        return coverage_dictionary

    def compute_size(self, map_dictionary):
        '''
        Takes a map_dictionary and returns a dictionary of sizes:
        {chrom: {polarity: {size: nbre of reads}}}
        '''
        size_dictionary = defaultdict(lambda: defaultdict(
                                      lambda: defaultdict(int)))
        #  to track empty chromosomes
        for chrom in self.chromosomes:
            if self.bam_object.count(chrom) == 0:
                size_dictionary[chrom]['F'][10] = 0
        for key in map_dictionary:
            for size in map_dictionary[key]:
                size_dictionary[key[0]][key[2]][size] += 1
        return size_dictionary

    def write_size_table(self, out):
        '''
        Dataset, Chromosome, Polarity, Size, Nbr_reads
        out is an *open* file handler
        '''
        for chrom in sorted(self.size):
            sizes = self.size[chrom]['F'].keys()
            sizes.extend(self.size[chrom]['R'].keys())
            for polarity in sorted(self.size[chrom]):
                for size in range(min(sizes), max(sizes)+1):
                    try:
                        line = [self.sample_name, chrom, polarity, size,
                                self.size[chrom][polarity][size]]
                    except KeyError:
                        line = [self.sample_name, chrom, polarity, size, 0]
                    line = [str(i) for i in line]
                    out.write('\t'.join(line) + '\n')

    def write_table(self, out):
        '''
        Dataset, Chromosome, Chrom_length, Coordinate, Nbr_reads
        Polarity, Max, Mean, Median, Coverage
        out is an *open* file handler
        '''
        for key in sorted(self.map_dict):
            line = [self.sample_name, key[0], self.chromosomes[key[0]],
                    key[1], len(self.map_dict[key]), key[2], self.max[key[0]],
                    self.mean[key], self.median[key], self.coverage[key]]
            line = [str(i) for i in line]
            out.write('\t'.join(line) + '\n')


def main(inputs, samples, file_out, size_file_out=''):
    F = open(file_out, 'w')
    header = ["Dataset", "Chromosome", "Chrom_length", "Coordinate",
              "Nbr_reads", "Polarity", "Max", "Mean", "Median", "Coverage"]
    F.write('\t'.join(header) + '\n')
    if size_file_out:
        Fs = open(size_file_out, 'w')
        header = ["Dataset", "Chromosome", "Polarity", "Size", "Nbr_reads"]
        Fs.write('\t'.join(header) + '\n')
        for file, sample in zip(inputs, samples):
            mapobj = Map(file, sample, computeSize=True)
            mapobj.write_table(F)
            mapobj.write_size_table(Fs)
        Fs.close()
    else:
        for file, sample in zip(inputs, samples):
            mapobj = Map(file, sample, computeSize=False)
            mapobj.write_table(F)
        F.close()


if __name__ == "__main__":
    args = Parser()
    # if identical sample names
    if len(set(args.sample_name)) != len(args.sample_name):
        args.sample_name = [name + '_' + str(i) for
                            i, name in enumerate(args.sample_name)]
    main(args.input, args.sample_name, args.output, args.sizes)
