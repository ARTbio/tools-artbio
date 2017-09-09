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

    def __init__(self, bam_file):
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.alignement_dic = self.index_alignments(self.bam_object)
        self.all_query_positions = self.query_positions(self.bam_object, overlap=10)
        self.readdic = self.make_readdic(self.bam_object)

    def make_readdic(self, bam_object):
        readdic = defaultdict(int)
        for read in bam_object.fetch():
            readdic[read.query_sequence] += 1
        return readdic

    def index_alignments(self, bam_object):
        '''
        dic[(chrom, pos, polarity)]: [readseq1, readseq2, ...] -> converted in set
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
            if pol == 'F' and len(self.alignement_dic[(chrom, coord+overlap-1, 'R')]) > 0:
                all_query_positions[chrom].append(coord)
        for chrom in all_query_positions:
            all_query_positions[chrom] = sorted(
                list(set(all_query_positions[chrom])))
        return all_query_positions
        
    def pairing(self, minquery, maxquery, mintarget, maxtarget, file,
                overlap=10):
        F = open(file, 'w')
        query_range = range(minquery, maxquery + 1)
        target_range = range(mintarget, maxtarget + 1)
        stringresult = []
        for chrom in self.all_query_positions:
            for pos in self.all_query_positions[chrom]:
               uppers = self.alignement_dic[chrom][pos]
               lowers = self.alignement_dic[chrom][pos+overlap-1]
               if uppers and lowers:
                   for upread in uppers:
                       for downread in lowers:
                           if (len(upread) in query_range and len(downread) in target_range) or (len(upread) in target_range and len(downread) in query_range):
                                stringresult.append(
                                    '>%s|%s|%s|%s|n=%s\n%s\n' %
                                    (chrom, pos+1, 'F', len(upread),
                                     self.readdic[upread], upread))
                                stringresult.append(
                                    '>%s|%s|%s|%s|n=%s\n%s\n' %
                                    (chrom, pos+overlap-len(downread), 'R', len(downread),
                                     self.readdic[downread], downread))
        stringresult = sorted(set(stringresult),
                              key=lambda x: stringresult.index(x))
        F.write(''.join(stringresult))

    def direct_pairing(self, minquery, maxquery, mintarget, maxtarget,
                       file, overlap=10):
        F = open(file, 'w')
        query_range = range(minquery, maxquery + 1)
        target_range = range(mintarget, maxtarget + 1)
        stringresult = []
        for chrom in sorted(self.chromosomes):
            for pos in (self.all_query_positions[chrom]):
                iterreads_1 = self.bam_object.fetch(chrom,
                                                    start=pos, end=pos+overlap-1)
                iterreads_2 = self.bam_object.fetch(chrom,
                                                    start=pos, end=pos+overlap-1)
                iterreads_3 = self.bam_object.fetch(chrom,
                                                    start=pos, end=pos+overlap-1)
                iterreads_4 = self.bam_object.fetch(chrom,
                                                    start=pos, end=pos+overlap-1)
                #  1
                for queryread in iterreads_1:
                    if queryread.reference_start == pos and \
                        queryread.query_alignment_length in query_range \
                            and not queryread.is_reverse:
                        for targetread in iterreads_2:
                            if (targetread.
                                get_reference_positions()[-1]
                                == queryread.get_reference_positions(
                                   )[overlap-1] and
                                targetread.query_alignment_length in
                                    target_range and targetread.is_reverse):
                                targetreadseq = self.revcomp(
                                    targetread.query_sequence)
                                stringresult.append(
                                    '>%s|%s|%s|%s|n=%s\n%s\n' %
                                    (chrom, queryread.reference_start+1,
                                     'F', queryread.query_alignment_length,
                                     self.readdic[queryread.query_sequence],
                                     queryread.query_sequence))
                                stringresult.append(
                                    '>%s|%s|%s|%s|n=%s\n%s\n' %
                                    (chrom, targetread.reference_start+1,
                                     'R', targetread.query_alignment_length,
                                     self.readdic[targetread.query_sequence],
                                     targetreadseq))
                #  2
                for queryread in iterreads_3:
                    if queryread.reference_end-1 == pos+overlap-1 and \
                        queryread.query_alignment_length in query_range \
                            and queryread.is_reverse:
                        for targetread in iterreads_4:
                            if (targetread.
                                reference_start
                                == pos and targetread.query_alignment_length
                                    in target_range and not
                                    targetread.is_reverse):
                                queryreadseq = self.revcomp(
                                    queryread.query_sequence)
                                targetreadseq = targetread.query_sequence
                                stringresult.append(
                                    '>%s|%s|%s|%s|n=%s\n%s\n' %
                                    (chrom, queryread.reference_start+1, 'R',
                                     queryread.query_alignment_length,
                                     self.readdic[queryread.query_sequence],
                                     queryreadseq))
                                stringresult.append(
                                    '>%s|%s|%s|%s|n=%s\n%s\n' %
                                    (chrom, targetread.reference_start+1,
                                     'F', targetread.query_alignment_length,
                                     self.readdic[targetread.query_sequence],
                                     targetreadseq))
        stringresult = sorted(set(stringresult),
                              key=lambda x: stringresult.index(x))
        F.write(''.join(stringresult))

    def revcomp(self, sequence):
        antidict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        revseq = sequence[::-1]
        return "".join([antidict[i] for i in revseq])


def main(input, minquery, maxquery, mintarget, maxtarget, output, overlap=10):
    mapobj = Map(input)
    mapobj.direct_pairing(minquery, maxquery, mintarget, maxtarget,
                          output, overlap)


if __name__ == "__main__":
    args = Parser()
    main(args.input, args.minquery, args.maxquery, args.mintarget,
         args.maxtarget, args.output, args.overlap)
