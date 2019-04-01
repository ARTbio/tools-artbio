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
    the_parser.add_argument('--bed', dest='bed', required=False,
                            help='Name of bed output must be specified\
                            if --cluster option used')
    the_parser.add_argument('--bed_skipsize', dest='bed_skipsize',
                            required=False, type=int, default=1,
                            help='Skip clusters of size equal or less than\
                            specified integer in the bed output. \
                            Default = 0, not skipping')
    the_parser.add_argument('--bed_skipdensity', dest='bed_skipdensity',
                            required=False, type=float, default=0,
                            help='Skip clusters of density equal or less than\
                            specified float number in the bed output. \
                            Default = 0, not skipping')
    the_parser.add_argument('--bed_skipcounts', dest='bed_skipcounts',
                            required=False, type=int, default=1,
                            help='Skip clusters of size equal or less than\
                            specified integer in the bed output. \
                            Default = 0, not skipping')
    the_parser.add_argument('--outputs', nargs='+', action='store',
                            help='list of two output paths (only two)')
    the_parser.add_argument('-M', '--plot_methods', nargs='+', action='store',
                            help='list of 2 plot methods (only two) among:\
                            Counts, Max, Mean, Median, Coverage and Size')
    the_parser.add_argument('--nostrand', action='store_true',
                            help='Consider reads regardless their polarity')

    args = the_parser.parse_args()
    return args


class Map:

    def __init__(self, bam_file, sample, minsize, maxsize, cluster, nostrand):
        self.sample_name = sample
        self.minsize = minsize
        self.maxsize = maxsize
        self.cluster = cluster
        if not nostrand:
            self.nostrand = False
        else:
            self.nostrand = True
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object, self.minsize,
                                        self.maxsize, self.nostrand)
        if self.cluster:
            self.map_dict = self.tile_map(self.map_dict, self.cluster)

    def create_map(self, bam_object, minsize, maxsize, nostrand=False):
        '''
        Returns a map_dictionary {(chromosome,read_position,polarity):
                                                    [read_length, ...]}
        '''
        map_dictionary = defaultdict(list)
        for chrom in self.chromosomes:
            # get empty value for start and end of each chromosome
            map_dictionary[(chrom, 1, 'F')] = []
            map_dictionary[(chrom, self.chromosomes[chrom], 'F')] = []
            if not nostrand:
                for read in bam_object.fetch(chrom):
                    positions = read.positions  # a list of covered positions
                    if read.is_reverse:
                        map_dictionary[(chrom, positions[-1]+1, 'R')].append(
                                        read.query_alignment_length)
                    else:
                        map_dictionary[(chrom, positions[0]+1, 'F')].append(
                                        read.query_alignment_length)
            else:
                for read in bam_object.fetch(chrom):
                    positions = read.positions  # a list of covered positions
                    if read.is_reverse:
                        map_dictionary[(chrom, positions[-1]+1, 'F')].append(
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
        and returns a map_dictionary with structure:
        {(chromosome,read_position,polarity):
            [*counts*, [start_clust, end_clust]]}
        '''
        clustered_dic = defaultdict(list)
        for chrom in self.chromosomes:
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
            # now 2 dictionnaries (F and R) with structure:
            # {centered_coordinate: [coord1, coord2, coord3, ..]}
            F_clust_dic = dict(zip(F_clust_keys, F_clust_values))
            R_clust_dic = dict(zip(R_clust_keys, R_clust_values))
            for centcoor in F_clust_dic:
                accumulator = []
                for coor in F_clust_dic[centcoor]:
                    accumulator.extend(map_dic[(chrom, coor, 'F')])
                '''
                compute the offset of the cluster due to
                size of reads
                '''
                last = sorted(F_clust_dic[centcoor])[-1]
                try:
                    margin = max(map_dic[(chrom, last, 'F')]) - 1
                except ValueError:
                    margin = 0
                clustered_dic[(chrom, centcoor, 'F')] = [len(accumulator), [
                    F_clust_dic[centcoor][0],
                    F_clust_dic[centcoor][-1] + margin]]
            for centcoor in R_clust_dic:
                accumulator = []
                for coor in R_clust_dic[centcoor]:
                    accumulator.extend(map_dic[(chrom, coor, 'R')])
                '''
                compute the offset of the cluster due to
                size of reads
                '''
                first = sorted(R_clust_dic[centcoor])[0]
                try:
                    margin = max(map_dic[(chrom, first, 'R')]) - 1
                except ValueError:
                    margin = 0
                clustered_dic[(chrom, centcoor, 'R')] = [len(accumulator), [
                    R_clust_dic[centcoor][0] - margin,
                    R_clust_dic[centcoor][-1]]]
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

    def compute_coverage(self, map_dictionary, out, quality=15):
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
                size_dictionary[chrom]['F'][10] = 0
        for key in map_dictionary:
            for size in map_dictionary[key]:
                size_dictionary[key[0]][key[2]][size] += 1
        self.write_size_table(size_dictionary, out)

    def write_table(self, mapdict, out):
        '''
        Writer of a tabular file
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
        Writer of a tabular file
        Dataset, Chromosome, Chrom_length, <category (size)>, <some value>
        out is an *open* file handler
        '''
        for chrom in sorted(sizedic):
            sizes = sizedic[chrom]['F'].keys()
            sizes.extend(sizedic[chrom]['R'].keys())
            for polarity in sorted(sizedic[chrom]):
                for size in range(min(sizes), max(sizes)+1):
                    try:
                        line = [self.sample_name, chrom, polarity, size,
                                sizedic[chrom][polarity][size]]
                    except KeyError:
                        line = [self.sample_name, chrom, polarity, size, 0]
                    line = [str(i) for i in line]
                    out.write('\t'.join(line) + '\n')

    def write_cluster_table(self, clustered_dic, out, bedpath):
        '''
        Writer of a tabular file
        Dataset, Chromosome, Chrom_length, Coordinate, Polarity,
        <some mapped value>
        out is an *open* file handler
        bed is an a file handler internal to the function
        '''
        def filterCluster(size, count, density):
            if size < args.bed_skipsize:
                return False
            if count < args.bed_skipcounts:
                return False
            if density <= args.bed_skipdensity:
                return False
            return True
        bed = open(bedpath, 'w')
        clusterid = 0
        for key in sorted(clustered_dic):
            start = clustered_dic[key][1][0]
            end = clustered_dic[key][1][1]
            size = end - start + 1
            read_count = clustered_dic[key][0]
            if self.nostrand:
                polarity = '.'
            elif key[2] == 'F':
                polarity = '+'
            else:
                polarity = '-'
            density = float(read_count) / size
            line = [self.sample_name, key[0], self.chromosomes[key[0]],
                    key[1], key[2], read_count,
                    str(start) + "-" + str(end), str(size), str(density)]
            line = [str(i) for i in line]
            out.write('\t'.join(line) + '\n')
            if filterCluster(size, read_count, density):
                clusterid += 1
                name = 'cluster_' + str(clusterid)
                bedline = [key[0], str(start-1), str(end), name,
                           str(read_count), polarity, str(density)]
                bed.write('\t'.join(bedline) + '\n')
        print("number of reported clusters:", clusterid)
        bed.close()


def main(inputs, samples, methods, outputs, minsize, maxsize, cluster,
         nostrand, bedfile=None, bed_skipsize=0):
    for method, output in zip(methods, outputs):
        out = open(output, 'w')
        if method == 'Size':
            header = ["Dataset", "Chromosome", "Polarity", method, "Counts"]
        elif cluster:
            header = ["Dataset", "Chromosome", "Chrom_length", "Coordinate",
                      "Polarity", method, "Start-End", "Cluster Size",
                      "density"]
        else:
            header = ["Dataset", "Chromosome", "Chrom_length", "Coordinate",
                      "Polarity", method]
        out.write('\t'.join(header) + '\n')
        for input, sample in zip(inputs, samples):
            mapobj = Map(input, sample, minsize, maxsize, cluster, nostrand)
            token = {"Counts": mapobj.compute_readcount,
                     "Max": mapobj.compute_max,
                     "Mean": mapobj.compute_mean,
                     "Median": mapobj.compute_median,
                     "Coverage": mapobj.compute_coverage,
                     "Size": mapobj.compute_size,
                     "cluster": mapobj.write_cluster_table}
            if cluster:
                token["cluster"](mapobj.map_dict, out, bedfile)
            else:
                token[method](mapobj.map_dict, out)
        out.close()


if __name__ == "__main__":
    args = Parser()
    # if identical sample names
    if len(set(args.sample_names)) != len(args.sample_names):
        args.sample_names = [name + '_' + str(i) for
                             i, name in enumerate(args.sample_names)]
    main(args.inputs, args.sample_names, args.plot_methods, args.outputs,
         args.minsize, args.maxsize, args.cluster, args.nostrand, args.bed)
