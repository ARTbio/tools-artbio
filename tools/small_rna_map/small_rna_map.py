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

	for read in bamfile:
		polarity = "Forward" if (not read.is_reverse) else "Reverse"
		if not read.is_unmapped:
			read_rname = read.reference_name
			read_position = read.pos + 1
			the_key = (read_rname, read_position)
			if ref_coord.has_key(the_key): 
				nmbr_read = ref_coord[the_key][0] + 1
				ref_coord[the_key]= (nmbr_read, polarity)		
			else :	
				ref_coord[the_key] = (1,polarity)
							
		continue
	
	return ref_coord #dictionary {(read_rname,read_position):(nmbr_read,polarity)}



def main(infile, output):
	with pysam.AlignmentFile(infile,"rb") as infile:
		with open(output,"a") as outfile:
			outfile.write("chromo" + "\t" + "coord" + "\t" + "nbr-reads" + "\t" + "polarity" + "\n")
			ref_coord = Ref_coord_pol(infile)
			ref_length = Ref_length(infile)
			for read in ref_coord:
				outfile.write(str(read[0]) + "\t" + str(ref_length[str(i[0])]) + "\t" + str(read[1]) + "\t" + str(ref_coord[read][0]) + "\t" + str(ref_coord[read][1])+ "\n")

if __name__ == "__main__":
	args= Parser()
	main(args.input, args.output)

