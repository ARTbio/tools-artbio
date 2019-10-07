import argparse
from collections import defaultdict

import pysam


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--inputs', dest='inputs', required=True,
                            nargs='+', help='list of input BAM files')
    the_parser.add_argument('--minsize', dest='minsize', type=int,
                            default=19, help='minimal size of reads')
    the_parser.add_argument('--maxsize', dest='maxsize', type=int,
                            default=29, help='maximal size of reads')
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
    the_parser.add_argument('--outputs', action='store',
                            help='list of two output paths (only two)')
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
        self.map_dict = self.create_map(self.bam_object, self.nostrand)
        if self.cluster:
            self.map_dict = self.tile_map(self.map_dict, self.cluster)

    def create_map(self, bam_object, nostrand=False):
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


def main(inputs, samples, outputs, minsize, maxsize, cluster,
         nostrand, bedfile=None, bed_skipsize=0):
    out = open(outputs, 'w')
    header = ["# Dataset", "Chromosome", "Chrom_length", "Coordinate",
              "Polarity", "Counts", "Start-End", "Cluster Size", "density"]
    out.write('\t'.join(header) + '\n')
    for input, sample in zip(inputs, samples):
        mapobj = Map(input, sample, minsize, maxsize, cluster, nostrand)
        mapobj.write_cluster_table(mapobj.map_dict, out, bedfile)
    out.close()


if __name__ == "__main__":
    args = Parser()
    # if identical sample names
    if len(set(args.sample_names)) != len(args.sample_names):
        args.sample_names = [name + '_' + str(i) for
                             i, name in enumerate(args.sample_names)]
    main(args.inputs, args.sample_names, args.outputs,
         args.minsize, args.maxsize, args.cluster, args.nostrand, args.bed)
