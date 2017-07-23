#! /bin/bash
# Execution command : bash small_rna_map.bash <sam file>
samtools view -h -o input.sam ../test-data/input.bam
samfile="input.sam"
cut -f2,3,4,5,6 ../test-data/output.tab > test-data.dat
grep -e "^@SQ" $samfile > header.tab
cat $samfile | perl -ne 'print unless /^@/ or ( (split /\t/)[1] != 0 and (split /\t/)[1] != 16 )' > output1.tab

cat << EOF > pyscript.py
lendic = {}
for line in open("header.tab"):
  fields = line[:-1].split()
  lendic[fields[1][3:]] = fields[2][3:]
readdic = {}
for line in (open("output1.tab")):
  fields = line[:-1].split()
  key = "-".join([fields[2],fields[3],"F" if fields[1]=="0" else "R"])
  try:
    readdic[key] += 1
  except KeyError:
    readdic[key] = 1
Out = open("output.tab", "w")
Out.write("Chromosome\tChrom_length\tCoordinate\tNbr_reads\tPolarity\n")
for key in sorted(readdic):
  fields = key.split("-")
  Chromosome, Coordinate, Polarity = fields[0], fields[1], fields[2]
  Chrom_length = lendic[Chromosome]
  Nbr_reads = str(readdic[key])
  Out.write("\t".join([Chromosome, Chrom_length, Coordinate, Nbr_reads, Polarity]) )
  Out.write("\n")
Out.close
EOF

python ./pyscript.py
rm output1.tab header.tab pyscript.py
sort -k 1,1 -k 3,3n -k 5,5  output.tab -o output.tab
if ! diff -q output.tab test-data.dat &>/dev/null; then
  >&2 echo "different"
else
  >&2 echo "Test passed"
fi

rm output.tab input.sam test-data.dat
