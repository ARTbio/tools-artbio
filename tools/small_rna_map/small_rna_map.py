import pysam
samfile= pysam.AlignmentFile("/home/artbio/Bureau/test.sam", "r")
output = open("/home/artbio/Bureau/output.gff", "w")
i=0
for read in samfile.fetch():
	if read.is_unmapped: NULL
	else: output.write(read.qname+"\t"+str(read.pos)+"\n")

output.close()

output 2= open("/home/artbio/Bureau/output2.gff", "w")
header_keys = samfile.header.keys()
Ref = samfile.references
lengths = samfile.lengths
i=0
for read in Ref:
	
	output2.write(read+"\t"+str(lengths[i])+"\n")
	i+=1

output2.close()
samfile.close()
