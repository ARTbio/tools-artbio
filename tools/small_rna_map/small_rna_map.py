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


def compute_map_dictionary(pysam_object):
    '''
    returns
    - map_dictionary {(read_chromosome,read_position,polarity):
                     [list_of_read_length]}
    - max_dictionary {chromosome_name:
                      max_of_number_of_read_at_any_position}
    '''
    unmatched_chromosomes = list(pysam_object.references)
    map_dictionary = defaultdict(list)
    max_dictionary = dict.fromkeys(pysam_object.references, 0)
    for read in pysam_object:
        if not read.is_unmapped:
            read_chromosome = read.reference_name
            # to shift from 0-based to 1-based coordinate
            read_position = read.pos + 1
            polarity = "F" if (not read.is_reverse) else "R"
            map_dictionary[(read_chromosome, read_position,
                            polarity)].append(read.qlen)
            if read_chromosome in unmatched_chromosomes:
                unmatched_chromosomes.remove(read_chromosome)
    for name in unmatched_chromosomes:
        # to put at least a null value by chromosome
        map_dictionary[(name, 1, 'F')] = [0]
    # make a dictionary with max value of reads by chromosome
    max_dictionary = defaultdict(int)
    for name in pysam_object.references:
        for key in map_dictionary.keys():
            if name in key and len(map_dictionary[key]) > max_dictionary[name]:
                max_dictionary[name] = len(map_dictionary[key])
    return map_dictionary, max_dictionary


class Map:

    def __init__(self, bam_file, sample):
        self.sample_name = sample
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = self.bam_object.references
        self.map_dict = self.create_map(self.bam_object)
        self.absent_chromosomes = self.get_absent_chrom(self.map_dict)
        self.max = self.compute_max(self.map_dict)

    def create_map(self, bam_object):
        '''
        Returns a map_dictionary {(read_chromosome,read_position,polarity):
                                                       [read_length, ...]}
        alternate method without empty chromosomes
        '''
        map_dictionary = defaultdict(list)
        for chrom in bam_object.references:
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
        for key in mean_dictionary:
            mean_dictionary = mean(mean_dictionary[key])
        return mean_dictionary


    def compute_median(self, map_dictionary):
        '''
        takes a map_dictionary as input and returns
        a mean_dictionary {(read_chromosome,read_position,polarity):
                                                       mean_value_of_reads}
        '''
        median_dictionary = dict()
        for key in median_dictionary:
            median_dictionary = median(median_dictionary[key])
        return median_dictionary


def main(inputs, sample_names, output):
    with open(output, 'w') as outfile:
        outfile.write('\t'.join(['Dataset', 'Chromosome', 'Chrom_length',
                                 'Coordinate', 'Nbr_reads', 'Polarity',
                                 'Max', 'Mean', 'Median\n']))
        chr_lengths = dict(zip(pysam.AlignmentFile(inputs[0], 'rb').references,
                               pysam.AlignmentFile(inputs[0], 'rb').lengths))
        for file, sample in zip(inputs, sample_names):
            with pysam.AlignmentFile(file, 'rb') as bamfile:
                map, maximum = compute_map_dictionary(bamfile)
                for record in sorted(map):
                    line = [sample, record[0],
                            str(chr_lengths[record[0]]),
                            str(record[1]), str(len(map[record])),
                            record[2], str(maximum[record[0]]),
                            str(round(numpy.mean(map[record]), 1)),
                            str(numpy.median(map[record]))]
                    outfile.write("\t".join(line) + "\n")


if __name__ == "__main__":
    args = Parser()
    # if identical sample names # to be tested
    if len(set(args.sample_name)) != len(args.sample_name):
        args.sample_name = [name + '_' + str(i) for
                            i, name in enumerate(args.sample_name)]
    main(args.input, args.sample_name, args.output)
