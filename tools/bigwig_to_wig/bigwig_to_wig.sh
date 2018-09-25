#!/bin/bash
#$Id: bigwig2wig 23 2014-01-28 12:09:22Z jens $

#SCRIPT CONVERTS BIGWIG FILE TO FIXED-STEP WIGGLE FORMAT FILE
#RESOLUTION IS CONTROLLED THROUGH THE BIN SIZE

#default bin_size
bin_size=500
mylab="wiggle file"

#parse input
while getopts hf:b:l: myarg
do	case "$myarg" in
	h)	echo "Usage: bigwig_correlation -f <bigwig_file> -b <bin_size>"
	        echo "Ex: bigwig_correlation -f <MYFILE.bw> -b 600"
		exit ;;
	f)	bigwig_file="$OPTARG" ;; #required
	l)	mylab="$OPTARG" ;; #optional 
	b)	bin_size="$OPTARG" ;; #optional
	[?])	echo "Usage: bigwig_correlation -f <MYFILE.bw> -b <bin_size>"
		exit 1 ;;
	esac
done

###################################################
###VALIDATE INPUT
###################################################

#make tmp-filename to hold chromosome info
org_assembly_file=$(mktemp -u)
bigWigInfo -chroms $bigwig_file | perl -ne "/^\tchr/ && print" | perl -pe "s/ +/\t/g" | cut -f2,4 > $org_assembly_file

#check bin_size & define step_size
bin_size_mod=$((bin_size % 2))  #determine modulus
if [ $bin_size_mod -ne 0 ]; then
  echo "Chosen bin_size must be an even positive number, added +1 to bin_size"
  bin_size=$((bin_size + 1))
fi

if [ $bin_size -lt 100 ]; then
  echo "ERROR: Chosen bin_size must be a positive number >=100"
  exit 1
fi
#set stetp size equal to bin size i.e. non-overlapping intervals
step_size=$bin_size

###################################################
###EXTRACT DENSITIES FROM NORMALIZED BIGWIG FILES
###################################################


#make track definition line
echo "track type=wiggle_0 name=$mylab description=\"fixedStep format\"" 

#for each chromsome
while read line; do
 
  cur_chr=$(echo $line | cut --delimiter=" " -f1)
  cur_length=$(echo $line | cut --delimiter=" " -f2)
  
  n_bins=$(echo "scale=0; (${cur_length}-${step_size})/${bin_size}" | bc)
 
  start=1
  stop=$(echo "$n_bins * $bin_size" | bc)

  #write header line for each chromosome
  echo "fixedStep chrom=$cur_chr start=$start step=$step_size span=$step_size"
  
  #get densities along chr in n_bins with chosen bin_size and step_size (giving overlap in bins)
  nice bigWigSummary $bigwig_file $cur_chr $start $stop $n_bins | perl -pe 's/\t/\n/g' | perl -pe "s/n\/a/0/" 
  #gives warning if no data in/for current chromosome
  
done < $org_assembly_file

#rm tmp
rm $org_assembly_file

