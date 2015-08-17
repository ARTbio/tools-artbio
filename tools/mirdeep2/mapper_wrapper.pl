#!/usr/bin/perl

$bowtie_index_name = $ARGV[0];
$output_arf = $ARGV[1];
$bowtie_files_path = $ARGV[2];
$options = join (" ", @ARGV[3..$#ARGV]);

$basename = `basename $bowtie_index_name`;
chomp $basename;

system ("mapper.pl $options -p $bowtie_files_path/$basename -t $output_arf 2> /dev/null");

if (-s $output_arf == 0) {die "No reads aligned to the reference.";}
