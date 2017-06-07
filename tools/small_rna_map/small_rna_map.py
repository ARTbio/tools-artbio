import os
from os.path import basename
import pysam
import argparse
import numpy
import collections


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--input', dest="input", required=True, nargs='+', help='input BAM files')
    the_parser.add_argument('--output', action="store", type=str, help="output tabular file")
    args = the_parser.parse_args()
    return args


def Ref_length(bamfile):
    """ 
    bamfile.references: list of reference names, and bamfile.lengths is list of their lengths 
    """
    ref_length={}
    for e1, e2 in zip(bamfile.references, bamfile.lengths):
        ref_length[e1]=e2
    return ref_length


def Ref_coord_pol(bamfile):
    ref_coord ={}
    for read in bamfile:
        if not read.is_unmapped:
            polarity = "F" if (not read.is_reverse) else "R"
            read_rname = read.reference_name
            read_position = read.pos + 1
            the_key = (read_rname, read_position,polarity)
            read_length = read.qlen
            try:
                ref_coord[the_key][0] = ref_coord[the_key][0]+1
                ref_coord[the_key][1].append(read_length)
            except:
                ref_coord[the_key] = [1,[read_length]]
    return collections.OrderedDict(sorted(ref_coord.items())) #dictionary {(read_rname,read_position,polarity):[nmbr_read,list_of_alignment_length]}


def calcul_stat(dictionary,ref):
    dico = {}	
    for keys in dictionary:
        key = keys[0] #reference name
        value = dictionary[keys][0] # position
        try:
            dico[key].append(value)
        except: 
            dico[key] = []
            dico[key].append(value)
    the_max = max(dico[ref])
    the_mean = numpy.mean(dico[ref])
    return (the_max, round(the_mean,2))


def filename(infile): 
    filename = basename(infile.filename)
    the_name = filename.split('.')
    return the_name[0]


def main(infile, output):
   with open(output,"w") as outfile:
       outfile.write("Dataset" + "\t" + "Chromosome"+ "\t" + "Chrom_length" + "\t" + "Coordinate" + "\t" + "Nbr_reads" + "\t" + "Polarity" + "\t" + "Max" + "\t" + "Mean" + "\t" + "Median" + "\n")
       for bamfile in infile:
            with pysam.AlignmentFile(bamfile,"rb") as bamfile:
                ref_coord = Ref_coord_pol(bamfile)
                ref_length = Ref_length(bamfile)
                tabulation = "\t"
                for read in ref_coord:
                    table = (str(filename(bamfile)), str(read[0]), str(ref_length[str(read[0])]), str(read[1]), str(ref_coord[read][0]), str(read[2]), str(calcul_stat(ref_coord,read[0])[0]), str(calcul_stat(ref_coord,read[0])[1]), str(numpy.median(ref_coord[read][1])))
                    outfile.write(tabulation.join(table) + "\n")


if __name__ == "__main__":
    args= Parser()
    main(args.input, args.output)
