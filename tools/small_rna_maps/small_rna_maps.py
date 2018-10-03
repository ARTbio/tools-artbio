import argparse
from collections import defaultdict

import numpy

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--inputs', dest='inputs', required=True,
                            nargs='+', help='list of input BAM files')
    the_parser.add_argument('--minsize', dest='minsize', type=int,
                            default=0, help='minimal size of reads')
    the_parser.add_argument('--maxsize', dest='maxsize', type=int,
                            default=10000, help='maximal size of reads')
    the_parser.add_argument('--cluster', dest='cluster', type=int,
                            default=0, help='clustering distance')
    the_parser.add_argument('--sample_names', dest='sample_names',
                            required=True, nargs='+',
                            help='list of sample names')
    the_parser.add_argument('--outputs', nargs='+', action='store',
                            help='list of two output paths (only two)')
    the_parser.add_argument('-M', '--plot_methods', nargs='+', action='store',
                            help='list of 2 plot methods (only two) among:\
                            Counts, Max, Mean, Median, Coverage and Size')
    args = the_parser.parse_args()
    return args


class Map:

    def __init__(self, bam_file, sample, minsize, maxsize, cluster):
        self.sample_name = sample
        self.minsize = minsize
        self.maxsize = maxsize
        self.cluster = cluster
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object, self.minsize,
                                        self.maxsize)
        if self.cluster:
            self.map_dict = self.tile_map(self.map_dict, self.cluster)

    def create_map(self, bam_object, minsize, maxsize):
        '''
        Returns a map_dictionary {(chromosome,read_position,polarity):
                                                    [read_length, ...]}
        '''
        map_dictionary = defaultdict(list)
        for chrom in self.chromosomes:
            # get empty value for start and end of each chromosome
            map_dictionary[(chrom, 1, 'F')] = []
            map_dictionary[(chrom, self.chromosomes[chrom], 'F')] = []
            for read in bam_object.fetch(chrom):
                positions = read.positions  # a list of covered positions
                if read.is_reverse:
                    map_dictionary[(chrom, positions[-1]+1, 'R')].append(
                                    read.query_alignment_length)
                else:
                    map_dictionary[(chrom, positions[0]+1, 'F')].append(
                                    read.query_alignment_length)
        return map_dictionary

    def grouper(self, iterable, clust_distance):
        prev = None
        group = []
        for item in iterable:
            if not prev or item - prev <= clust_distance:
                group.append(item)
            else:
                yield group
                group = [item]
            prev = item
        if group:
            yield group

    def tile_map(self, map_dic, clust_distance):
        '''
        takes a map_dictionary {(chromosome,read_position,polarity):
                                                    [read_length, ...]}
        and retur a map_dictionary with same structure but with
        read positions aggregated by size
        '''
        clustered_dic = defaultdict(list)
        for chrom in self.chromosomes:
            clustered_dic[(chrom, 1, 'F')] = []
            clustered_dic[(chrom, self.chromosomes[chrom], 'F')] = []
            F_chrom_coord = []
            R_chrom_coord = []
            for key in map_dic:
                if key[0] == chrom and key[2] == 'F':
                    F_chrom_coord.append(key[1])
                elif key[0] == chrom and key[2] == 'R':
                    R_chrom_coord.append(key[1])
            F_chrom_coord = list(set(F_chrom_coord))
            R_chrom_coord = list(set(R_chrom_coord))
            F_chrom_coord.sort()
            R_chrom_coord.sort()
            F_clust_values = [i for i in self.grouper(F_chrom_coord,
                                                      clust_distance)]
            F_clust_keys = [(i[-1]+i[0])/2 for i in F_clust_values]
            R_clust_values = [i for i in self.grouper(R_chrom_coord,
                                                      clust_distance)]
            R_clust_keys = [(i[-1]+i[0])/2 for i in R_clust_values]
            F_clust_dic = dict(zip(F_clust_keys, F_clust_values))
            R_clust_dic = dict(zip(R_clust_keys, R_clust_values))
            # {centered_coordinate: [coord1, coord2, coord3, ..]}
            for centcoor in F_clust_dic:
                accumulator = []
                for coor in F_clust_dic[centcoor]:
                    accumulator.extend(map_dic[(chrom, coor, 'F')])
                clustered_dic[(chrom, centcoor, 'F')] = accumulator
            for centcoor in R_clust_dic:
                accumulator = []
                for coor in R_clust_dic[centcoor]:
                    accumulator.extend(map_dic[(chrom, coor, 'R')])
                clustered_dic[(chrom, centcoor, 'R')] = accumulator
        return clustered_dic

    def compute_readcount(self, map_dictionary, out):
        '''
        takes a map_dictionary as input and writes
        a readmap_dictionary {(chromosome,read_position,polarity):
                              number_of_reads}
        in an open file handler out
        '''
        readmap_dictionary = dict()
        for key in map_dictionary:
            readmap_dictionary[key] = len(map_dictionary[key])
        self.write_table(readmap_dictionary, out)

    def compute_max(self, map_dictionary, out):
        '''
        takes a map_dictionary as input and writes
        a max_dictionary {(chromosome,read_position,polarity):
                              max_of_number_of_read_at_any_position}
        Not clear this function is still required
        '''
        merge_keylist = [(i[0], 0) for i in map_dictionary.keys()]
        max_dictionary = dict(merge_keylist)
        for key in map_dictionary:
            if len(map_dictionary[key]) > max_dictionary[key[0]]:
                max_dictionary[key[0]] = len(map_dictionary[key])
        self.write_table(max_dictionary, out)

    def compute_mean(self, map_dictionary, out):
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
        self.write_table(mean_dictionary, out)

    def compute_median(self, map_dictionary, out):
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
        self.write_table(median_dictionary, out)

    def compute_coverage(self, map_dictionary, out, quality=10):
        '''
        takes a map_dictionary as input and returns
        a coverage_dictionary {(chromosome,read_position,polarity):
                                                coverage}
        '''
        coverage_dictionary = dict()
        for chrom in self.chromosomes:
            coverage_dictionary[(chrom, 1, 'F')] = 0
            coverage_dictionary[(chrom, self.chromosomes[chrom], 'F')] = 0
            for read in self.bam_object.fetch(chrom):
                positions = read.positions  # a list of covered positions
                for pos in positions:
                    if not map_dictionary[(chrom, pos+1, 'F')]:
                        map_dictionary[(chrom, pos+1, 'F')] = []
        for key in map_dictionary:
            coverage = self.bam_object.count_coverage(
                                                reference=key[0],
                                                start=key[1]-1,
                                                end=key[1],
                                                quality_threshold=quality)
            """ Add the 4 coverage values """
            coverage = [sum(x) for x in zip(*coverage)]
            coverage_dictionary[key] = coverage[0]
        self.write_table(coverage_dictionary, out)

    def compute_size(self, map_dictionary, out):
        '''
        Takes a map_dictionary and returns a dictionary of sizes:
        {chrom: {polarity: {size: nbre of reads}}}
        '''
        size_dictionary = defaultdict(lambda: defaultdict(
                                      lambda: defaultdict(int)))
        #  to track empty chromosomes
        for chrom in self.chromosomes:
            if self.bam_object.count(chrom) == 0:
                size_dictionary[chrom] = 0
        for key in map_dictionary:
            for size in map_dictionary[key]:
                size_dictionary[key[0]][key[2]][size] += 1
        self.write_size_table(size_dictionary, out)

    def write_table(self, mapdict, out):
        '''
        Generic writer
        Dataset, Chromosome, Chrom_length, Coordinate, Polarity,
        <some mapped value>
        out is an *open* file handler
        '''
        for key in sorted(mapdict):
            line = [self.sample_name, key[0], self.chromosomes[key[0]],
                    key[1], key[2], mapdict[key]]
            line = [str(i) for i in line]
            out.write('\t'.join(line) + '\n')

    def write_size_table(self, sizedic, out):
        '''
        Generic writer of summary values
        Dataset, Chromosome, Chrom_length, <some category>, <some value>
        out is an *open* file handler
        '''
        for chrom in sorted(sizedic):
            try:
                sizes = sizedic[chrom]['F'].keys()
                sizes.extend(sizedic[chrom]['R'].keys())
                for polarity in sorted(sizedic[chrom]):
                    for size in range(min(sizes), max(sizes)+1):
                        line = [self.sample_name, chrom, polarity, size,
                                sizedic[chrom][polarity][size]]
                        line = [str(i) for i in line]
                        out.write('\t'.join(line) + '\n')
            except TypeError:
                # Print unmapped chromosomes so as to get empty graphs
                line = [self.sample_name, chrom, 'F', 0, 0]
                line = [str(i) for i in line]
                out.write('\t'.join(line) + '\n')
                line = [self.sample_name, chrom, 'R', 0, 0]
                line = [str(i) for i in line]
                out.write('\t'.join(line) + '\n')


def main(inputs, samples, methods, outputs, minsize, maxsize, cluster):
    for method, output in zip(methods, outputs):
        F = open(output, 'w')
        if method == 'Size':
            header = ["Dataset", "Chromosome", "Polarity", method, "Counts"]
        else:
            header = ["Dataset", "Chromosome", "Chrom_length", "Coordinate",
                      "Polarity", method]
        F.write('\t'.join(header) + '\n')
        for input, sample in zip(inputs, samples):
            mapobj = Map(input, sample, minsize, maxsize, cluster)
            token = {"Counts": mapobj.compute_readcount,
                     "Max": mapobj.compute_max,
                     "Mean": mapobj.compute_mean,
                     "Median": mapobj.compute_median,
                     "Coverage": mapobj.compute_coverage,
                     "Size": mapobj.compute_size}
            token[method](mapobj.map_dict, F)
        F.close()


if __name__ == "__main__":
    args = Parser()
    # if identical sample names
    if len(set(args.sample_names)) != len(args.sample_names):
        args.sample_names = [name + '_' + str(i) for
                             i, name in enumerate(args.sample_names)]
    main(args.inputs, args.sample_names, args.plot_methods, args.outputs,
         args.minsize, args.maxsize, args.cluster)
