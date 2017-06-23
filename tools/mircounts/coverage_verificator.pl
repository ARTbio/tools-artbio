#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
##
# Takes alignment file in sam format and lattice dataframe and checks if the coverage is correct
# if it isn't it will print the chromosome and the position were the coverage isn't ok 
# the last column is the number of surplus reads if it's positif and missing reads if it's negatif
##
my $aln = '';
my $lattice = '';
my %hash = ();
my %lattice_hash = ();
GetOptions("alignment=s" => \$aln,
           "lattice=s"   => \$lattice)or die ("Error in command line arguments\n");

if(!$aln || !$lattice) {
  die ("Need all options alignment and lattice\n");
}


open(my $LAT, '<', $lattice)or die("Cannot open $lattice\n");
my $header =<$LAT>;
while(<$LAT>) {
  if(/^\S+\t(\S+)\t(\d+)\t\S+\t(\d+)/) {
    my $chrom = $1;
    my $offset = $2;
    my $count = $3;
    $lattice_hash{$chrom}{$offset} = $count;
  }
}
close($LAT);

open(my $AL, '<', $aln)or die("Cannot open $aln\n");
while(<$AL>) {
  if(/^[^@]/){
    chomp $_;
    my @line = split(/\t/,$_);
    my $chrom = $line[2];
    my $start = $line[3];
    my $stop = length($line[9]);
    for(my $i=$start;$i<$stop+1;$i++){
      if(!$lattice_hash{$chrom}{$i}){
        $i = $stop+1;
      }else{
        if ($chrom eq 'dme-miR-986-5p'){print 'it is '.$chrom."\tpos $i\tstop $stop\tnbreadsbeafor- $lattice_hash{$chrom}{$i}\n";}
        $lattice_hash{$chrom}{$i} -= 1;
      }
    }
  }
}
close($AL);

foreach my $chrom(keys(%lattice_hash)){
  foreach my $pos(sort(keys($lattice_hash{$chrom}))){
    if($lattice_hash{$chrom}{$pos} != 0) {
      print "$chrom\t$pos\t$lattice_hash{$chrom}{$pos}\n";
    }
  }
}
