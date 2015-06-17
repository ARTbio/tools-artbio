#!/bin/bash

min_k=$1
max_k=$2
input=$3
job_dir=$4

for i in `seq $min_k 2 $max_k`
do
     velveth $job_dir/outputFolder_$i $i $input
     velvetg $job_dir/outputFolder_$i -read_trkg yes
     oases $job_dir/outputFolder_$i
done

velveth $job_dir/MergedAssemblyFolder 27 -long outputFolder_*/transcripts.fa
velvetg $job_dir/MergedAssemblyFolder -read_trkg yes -conserveLong yes
oases $job_dir/MergedAssemblyFolder -merge yes

