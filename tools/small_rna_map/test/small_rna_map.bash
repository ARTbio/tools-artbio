#! /bin/bash
# Executing this script require bash v4
# Execution command : bash small_rna_map.bash

sort -k 3 -k 4  samfile.sam -o samfile.sam

while read Line
do 
   if [[ ${Line:0:1} != "@" ]] ; then
       echo "$Line" >> output1.tab
   elif [[ ${Line:0:3} = "@SQ" ]] ; then
           echo "$Line " >> header.tab
   fi
done < samfile.sam

while read qname flag rname pos mapq cigar mrnm mpos isize seq qual opt
do 
   if [[ $flag = 0 ]] ; then
       polarity="F"
   elif [[ $flag = 16 ]] ;
   then
       polarity="R"
   fi
   if [[ $flag = 0 ]] || [[ $flag = 16 ]] ; then
       echo  "$rname   $pos    $polarity   ${#seq}" >> output2.tab
   fi
done < output1.tab
rm output1.tab

while read SQ SN LN
do 
    declare -A tab
    tab["${SN:3}"]=${LN:3};
done < header.tab
rm header.tab
#echo "Nombre d'éléments : "${#tab[*]} ;
# rm header.tab

while read rname coord polarity qlen
do  
    key=$rname$coord$polarity
    len="${tab[$rname]}"
    echo "$rname    $len    $coord  $qlen   $polarity   $key" >> output3.tab
   
done < output2.tab
rm output2.tab

echo -e "Chromosome\tChrom_length\tCoordinate\tNbr_reads\tPolarity" > output.tab 

while read rname rlen coord qlen pol keys
do  

    if [[ $key != $keys ]] ; then
        key=$keys
        nbr_reads=`cat output3.tab | grep "$key" | wc -l`
        echo -e "$rname\t$rlen\t$coord\t$nbr_reads\t$pol" >> output.tab 
    fi
done < output3.tab
rm output3.tab

sort -k 1 -k 3 -k 4  output.tab -o output.tab

