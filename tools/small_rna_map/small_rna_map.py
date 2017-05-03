import pysam
import argparse
import numpy

def Parser():
        the_parser = argparse.ArgumentParser()
        the_parser.add_argument('--input', action="store", type=str, help="input BAM file")
        the_parser.add_argument('--output', action="store", type=str, help="output tabular file")

        args = the_parser.parse_args()
        return args

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

def calcul_stat(dictionary,ref):
	dico = {}
	
	for key,value in dictionary:
		if dico.has_key(key):
			dico[key].append(value)

		else: 
			dico[key] = []
			dico[key].append(value)



	the_max = max(dico[ref])
	the_mean = numpy.mean(dico[ref])
	the_median = numpy.median(dico[ref])
	return (the_max, the_mean, the_median)

def main(infile, output):
        with pysam.AlignmentFile(infile,"rb") as infile:
                with open(output,"a") as outfile:
                        outfile.write("chromo" + "\t" + "coord" + "\t" + "nbr-reads" + "\t" + "polarity" + "\t" + "Max" + "\t" + "Mean" + "\t" + "Median" + "\n")
                        ref_coord = Ref_coord_pol(infile)
                        ref_length = Ref_length(infile)
                        for read in ref_coord:
                                outfile.write(str(read[0]) + "\t" + str(ref_length[str(read[0])]) + "\t" + str(read[1]) + "\t" + str(ref_coord[read][0])+ "\t" + str(ref_coord[read][1]) + "\t" + str(calcul_stat(ref_coord,read[0])[0]) + "\t" + str(calcul_stat(ref_coord,read[0])[1]) + "\t" + str(calcul_stat(ref_coord,read[0])[2])+ "\n")

if __name__ == "__main__":
        args= Parser()
        main(args.input, args.output)


