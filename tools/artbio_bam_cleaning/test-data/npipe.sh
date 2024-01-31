input="chr22_sample" && \
genome="chr22.fa"
samtools index $input".bam" && \
sambamba view -h -t 8 --filter='mapping_quality >= 1 and not(unmapped) and not(mate_is_unmapped) and not(duplicate)' -f 'bam' ${input}".bam" \
| bamleftalign --fasta-reference $genome -c --max-iterations "5" - \
| samtools calmd  -C 50 -b -@ 8 - $genome > $input".filt1.dedup.bamleft.calmd.bam" && \
sambamba view -h -t 8 --filter='mapping_quality <= 254' -f 'bam' -o $input".filt1.dedup.bamleft.calmd.filt2.bam" $input".filt1.dedup.bamleft.calmd.bam"

