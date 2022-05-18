input="match_chr21_DBA_974" && \
genome="chr21.fa" && \
samtools index $input".bam" && \
sambamba view -h -t 8 --filter='mapping_quality >= 1 and not(unmapped) and not(mate_is_unmapped)' -f 'bam' $input".bam" \
| samtools rmdup - - \
tee $input".filt1.dedup.bam"| bamleftalign --fasta-reference $genome -c --max-iterations "5" - \
| samtools calmd  -C 50 -b -@ 4 - $genome > $input".filt1.dedup.bamleft.calmd.bam" && \
sambamba view -h -t 8 --filter='mapping_quality <= 254' -f 'bam' -o $input".filt1.dedup.bamleft.calmd.filt2.bam" $input".filt1.dedup.bamleft.calmd.bam"

