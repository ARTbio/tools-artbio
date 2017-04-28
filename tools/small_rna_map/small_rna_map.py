import pysam
import argparse
from numpy import mean, median

def Parser():
	the_parser = argparse.ArgumentParser()
	the_parser.add_argument('--input', action="store", type=str, help="input BAM file")
	the_parser.add_argument('--output', action="store", type=str, help="output tabular file")

	args = the_parser.parse_args()
	return args

bamfile= pysam.AlignmentFile("/home/artbio/Bureau/input.bam", "rb")
output = open("/home/artbio/Bureau/output2.tab", "w")

def Ref_length(bamfile):
	Ref = bamfile.references #Reference names
	lengths = bamfile.lengths # Reference lengths
	filename = bamfile.filename # filename
	ref_length = {}
	j=0
	for ref in Ref:
		ref_length[ref]=lengths[j]
		j += 1
	return ref_length		



def Ref_coord_pol(bamfile):
	ref_coord ={}
	i=0
	for read in bamfile:
		polarity = "Forward" if (not read.is_reverse) else "Reverse"
		if not read.is_unmapped:
			read_rname = read.reference_name
			read_position = read.pos
			
			if ref_coord.has_key(read_rname): 
				
				if int(ref_coord[read_rname][i])==read_position:
					ref_coord[read_rname][i] += 1 
				else : 
					ref_coord[read_rname].append(1)
					i = i+1				
			else :	
				ref_coord[read_rname] = [1]
				i=0			
		continue
		ref_coord[read_rname].append(polarity)
	return ref_coord



			
output.close()
bamfile.close()

def main(infile, output):
	
	print "\n"

if __name__ == "__main__":
	args= Parser()
	main(args.input, args.output)

