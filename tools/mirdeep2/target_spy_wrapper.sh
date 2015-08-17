export PATH=/share/apps/ViennaRNA-1.6.1/bin:/share/apps/TargetSpy/bin:$PATH

TargetSpy -microRNAs $1 -transcripts $2 -result targets
gunzip targets.gz
