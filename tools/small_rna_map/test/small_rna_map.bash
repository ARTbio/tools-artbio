#! /bin/bash

sort -k 3 -k 4  samfile.sam -o samfile.sam

while read Line
do 
   if [ ${Line:0:1} != "@" ] ; then
       echo "$Line" >> output1.tab
   elif [ ${Line:0:3} = "@SQ" ] ; then
           echo "$Line " >> header.tab
   fi
done < samfile.sam

while read qname flag rname pos mapq cigar mrnm mpos isize seq qual opt
do 
   if [ $flag = 0 ] ; then
       polarity="F"
   elif [ $flag = 16 ] ;
   then
       polarity="R"
   fi
   if [ $flag = 0 ] || [ $flag = 16 ] ; then
       echo  "$rname   $pos    $polarity   ${#seq}" >> output2.tab
   fi
done < output1.tab

while read SQ SN LN
do 
    declare -A tab
    tab["${SN:3}"]=${LN:3};
done < header.tab

#echo "Nombre d'éléments : "${#tab[*]} ;
# rm header.tab

while read rname coord polarity qlen
do  
    len="${tab[$rname]}"
    echo "$rname    $len    $coord  $qlen   $polarity" >> output3.tab
   
done < output2.tab

declare -A query
while read rname rlen coord qlen pol
do
    key=$rname$coord$pol
    query[$key]="${query[$key]}",$qlen
    echo "$key  ${query[$key]}  $rname  $coord  $pol" >> output4.tab


done < output3.tab


# cat output4.tab | grep "$key" | wc -l

# $(($nbr_read[$key]+1))
