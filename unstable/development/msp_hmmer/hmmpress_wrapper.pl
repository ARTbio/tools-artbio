#!/usr/bin/perl

use warnings;
use strict;

die("Infile, outfile, outdir required\n") unless @ARGV == 3;
my ($infile, $outfile, $outdir)=@ARGV;

# setup
-d $outdir or mkdir($outdir) or die($!);
symlink($infile,"$outdir/hmm");

# create output header file
open(IN, "<$infile") or die($!);
open(OUT, ">$outfile") or die($!);
my $ok=0;
while (my $line=<IN>) {
    if (!$ok) {
        die("Invalid input file (HMMER3 format required)\n") unless $line =~ /^HMMER3/;
        $ok=1;
    } elsif ($line =~ /^HMM/) {
        last;
    } else {
        print OUT $line;
    }
}
close IN;
close OUT;

# hmmpress
my $output=`hmmpress $outdir/hmm 2>&1`;
if ($? != 0) {
    $output="FAILED\n" unless $output;
    die($output);
}
unlink("$outdir/hmm");
my @output=split(/\n/, $output);
print $output[1], "\n";
exit;
