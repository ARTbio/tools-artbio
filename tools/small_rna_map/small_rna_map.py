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
