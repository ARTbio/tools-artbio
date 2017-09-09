import argparse
from collections import defaultdict

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
        '--overlap', type=int,
        help="Overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--output', action="store", type=str,
        help="Pairable sequences")
    args = the_parser.parse_args()
    return args


class Map:

    def __init__(self, bam_file, output, minquery=23, maxquery=29,
                 mintarget=23, maxtarget=29, overlap=10):
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.output = output
        self.query_range = range(minquery, maxquery + 1)
        self.target_range = range(mintarget, maxtarget + 1)
        self.overlap = overlap
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.alignement_dic = self.index_alignments(self.bam_object)
        self.all_query_positions = self.query_positions(self.bam_object,
                                                        overlap=self.overlap)
        self.readdic = self.make_readdic(self.bam_object)
        self.pairing()

    def make_readdic(self, bam_object):
        readdic = defaultdict(int)
        for read in bam_object.fetch():
            readdic[read.query_sequence] += 1
        return readdic

    def index_alignments(self, bam_object):
        '''
        dic[(chrom, pos, polarity)]: [readseq1, readseq2, ...]
        the list value is further converted in set
        '''
        dic = defaultdict(list)
        for chrom in self.chromosomes:
            for read in bam_object.fetch(chrom):
                if read.is_reverse:
                    coord = read.reference_end-1
                    pol = 'R'
                else:
                    coord = read.reference_start
                    pol = 'F'
                dic[(chrom, coord, pol)].append(read.query_sequence)
        for key in dic:
            dic[key] = set(dic[key])
        return dic

    def query_positions(self, bam_object, overlap):
        all_query_positions = defaultdict(list)
        for genomicKey in self.alignement_dic.keys():
            chrom, coord, pol = genomicKey
            if pol == 'F' and len(self.alignement_dic[(chrom,
                                                      coord+overlap-1,
                                                      'R')]) > 0:
                all_query_positions[chrom].append(coord)
        for chrom in all_query_positions:
            all_query_positions[chrom] = sorted(
                list(set(all_query_positions[chrom])))
        return all_query_positions

    def pairing(self):
        F = open(self.output, 'w')
        query_range = self.query_range
        target_range = self.target_range
        overlap = self.overlap
        stringresult = []
        header_template = '>%s|coord=%s|strand %s|size=%s|nreads=%s\n%s\n'
        for chrom in sorted(self.chromosomes):
            for pos in self.all_query_positions[chrom]:
                stringbuffer = []
                uppers = self.alignement_dic[chrom, pos, 'F']
                lowers = self.alignement_dic[chrom, pos+overlap-1, 'R']
                if uppers and lowers:
                    for upread in uppers:
                        for downread in lowers:
                            if (len(upread) in query_range and len(downread) in
                                target_range) or (len(upread) in target_range
                                                  and len(downread) in
                                                  query_range):
                                stringbuffer.append(
                                    header_template %
                                    (chrom, pos+1, '+', len(upread),
                                     self.readdic[upread], upread))
                                stringbuffer.append(
                                    header_template %
                                    (chrom, pos+overlap-len(downread)+1, '-',
                                     len(downread), self.readdic[downread],
                                     self.revcomp(downread)))
                stringresult.extend(sorted(set(stringbuffer)))
        F.write(''.join(stringresult))

    def revcomp(self, sequence):
        antidict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        revseq = sequence[::-1]
        return "".join([antidict[i] for i in revseq])


if __name__ == "__main__":
    args = Parser()
    mapobj = Map(args.input, args.output, args.minquery, args.maxquery,
                 args.mintarget, args.maxtarget, args.overlap)
