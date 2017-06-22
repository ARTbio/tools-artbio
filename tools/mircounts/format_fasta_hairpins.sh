GENOME_KEY=$1
        
gunzip hairpin.fa.gz
sed -i.bak '/^[^>]/ y/uU/tT/' hairpin.fa ## replace U by tT
sed -i.bak2 -E 's/ .+//' hairpin.fa ## just leaves mir name as one word header
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < hairpin.fa > hairpin.fa.bak3
tail -n +2 hairpin.fa.bak3 > hairpin.fa ## generate single line sequences
awk 'BEGIN{RS=">"}{gsub("\n","\t",$0); print ">"$0}' < hairpin.fa > hairpin.fa.tmp
mv hairpin.fa hairpin.bak4 && tail -n +2 hairpin.fa.tmp > hairpin.fa
rm hairpin.fa.tmp ## tabular sequences
sed -i.bak5 -E $'s/\t$//g' hairpin.fa ## remove tab before end line leaved by previous awk
grep ">${GENOME_KEY}-" hairpin.fa > hairpin.fa.tmp
mv hairpin.fa hairpin.fa.bak6
mv hairpin.fa.tmp hairpin.fa ## filter tabular hairpins with proper genomeKey
tr '\t' '\n' < hairpin.fa > hairpin.fa.tmp
mv hairpin.fa hairpin.fa.bak7
mv hairpin.fa.tmp hairpin.fa ## terminate parsing by regenerating fasta format, bowtie-build ready
rm ./*.bak* ## cleaning job directory

