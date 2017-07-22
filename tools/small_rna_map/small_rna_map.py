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
    the_parser.add_argument('--output', action="store",
                            type=str, help='output tabular file')
    args = the_parser.parse_args()
    return args


class Map:

    def __init__(self, bam_file, sample):
        self.sample_name = sample
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object)
        self.max = self.compute_max(self.map_dict)
        self.mean = self.compute_mean(self.map_dict)
        self.median = self.compute_median(self.map_dict)
        self.coverage = self.compute_coverage(self.map_dict) 

    def create_map(self, bam_object):
        '''
        Returns a map_dictionary {(read_chromosome,read_position,polarity):
                                                       [read_length, ...]}
        alternate method without empty chromosomes
        '''
        map_dictionary = dict()
        for chrom in self.chromosomes:
            for pos in range (self.chromosomes[chrom]):
                map_dictionary[(chrom, pos+1, 'F')] = []  # get all chromosomes
                map_dictionary[(chrom, pos+1, 'R')] = []  # get all chromosomes
        for chrom in self.chromosomes:
            for read in bam_object.fetch(chrom):
                read_positions = read.positions  # a list of covered positions
                if read.is_reverse:
                    map_dictionary[(chrom, read_positions[-1]+1,
                                    'R')].append(read.query_alignment_length)
                else:
                    map_dictionary[(chrom, read_positions[0]+1,
                                    'F')].append(read.query_alignment_length)
        return map_dictionary

    def compute_max(self, map_dictionary):
        '''
        takes a map_dictionary as input and returns
        a max_dictionary {chromosome_name:
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
        a mean_dictionary {(read_chromosome,read_position,polarity):
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
        a mean_dictionary {(read_chromosome,read_position,polarity):
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
        """
        Takes a AlignmentFile object and returns a dictionary of lists
        of coverage along the coordinates of pre_mirs (as keys)
        """
        coverage_dictionary = dict()
        for chrom in self.chromosomes:
            coverage = self.bam_object.count_coverage(reference=chrom,
                                                start=0,
                                                end=self.chromosomes[chrom],
                                                quality_threshold=quality)
            """ Add the 4 coverage values """
            coverage = [sum(x) for x in zip(*coverage)]
            for pos, cov in enumerate(coverage):
                coverage_dictionary[(chrom, pos+1, "F")] = cov
                coverage_dictionary[(chrom, pos+1, "R")] = cov
        return coverage_dictionary

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


def main(inputs, samples, file_out):
    F = open(file_out, 'w')
    header = ["Dataset", "Chromosome", "Chrom_length", "Coordinate",
              "Nbr_reads", "Polarity", "Max", "Mean", "Median", "Coverage"]
    header = '\t'.join(header) + '\n'
    F.write(header)
    for file, sample in zip(inputs, samples):
        mapobj = Map(file, sample)
        mapobj.write_table(F)
    F.close()


if __name__ == "__main__":
    args = Parser()
    # if identical sample names # to be tested
    if len(set(args.sample_name)) != len(args.sample_name):
        args.sample_name = [name + '_' + str(i) for
                            i, name in enumerate(args.sample_name)]
    main(args.input, args.sample_name, args.output)
