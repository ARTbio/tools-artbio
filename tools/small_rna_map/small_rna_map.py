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

Ref = bamfile.references #Reference names
lengths = bamfile.lengths # Reference lengths
nbr_read=0

for ref in Ref:
	for read in bamfile:
		if not read.is_unmapped:
			if ref==bamfile.getrname(read.reference_id):
				nbr_read=1+nbr_read
	
			
	output.write(ref+"\t"+str(read.rlen)+"\t"+str(nbr_read)+"\n")
			
output.close()
bamfile.close()

def main(infile, output):
	
	print "\n"

if __name__ == "__main__":
	args= Parser()
	main(args.input, args.output)

