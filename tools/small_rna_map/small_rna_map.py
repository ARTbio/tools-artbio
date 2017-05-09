import os
from os.path import basename
import pysam
import argparse
import numpy
import collections


def Parser():
        the_parser = argparse.ArgumentParser()
        the_parser.add_argument('--input', action="store", type=str, help="input BAM file")
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
                        try:
                                ref_coord[the_key] = ref_coord[the_key] + 1
                                
                        except:
                                ref_coord[the_key] = 1
        return collections.OrderedDict(sorted(ref_coord.items())) #dictionary {(read_rname,read_position,polarity):nmbr_read}


def calcul_stat(dictionary,ref):
	dico = {}	
	for keys in dictionary:
		key = keys[0] #reference name
		value = keys[1] # position
		try:
			dico[key].append(value)
		except: 
			dico[key] = []
			dico[key].append(value)
	the_max = max(dico[ref])
	the_mean = numpy.mean(dico[ref])
	the_median = numpy.median(dico[ref])
	return (the_max, round(the_mean,2), the_median)


def filename(infile): 
	filename = basename(infile.filename)
	the_name = filename.split('.')
	return the_name[0]


def main(infile, output):
        with pysam.AlignmentFile(infile,"rb") as infile:
                with open(output,"w") as outfile:
                        outfile.write("Dataset" + "\t" + "Chromosome"+ "\t" + "Chrom_length" + "\t" + "Coordinate" + "\t" + "Nbr_reads" + "\t" + "Polarity" + "\t" + "Max" + "\t" + "Mean" + "\t" + "Median" + "\n")
                        ref_coord = Ref_coord_pol(infile)
                        ref_length = Ref_length(infile)
                        for read in ref_coord:
				outfile.write(str(filename(infile)) + "\t" + str(read[0]) + "\t" + str(ref_length[str(read[0])]) + "\t" + str(read[1]) + "\t" + str(ref_coord[read])+ "\t" + str(read[2]) + "\t" + str(calcul_stat(ref_coord,read[0])[0]) + "\t" + str(calcul_stat(ref_coord,read[0])[1]) + "\t" + str(calcul_stat(ref_coord,read[0])[2])+ "\n")


if __name__ == "__main__":

        args= Parser()
        main(args.input, args.output)

