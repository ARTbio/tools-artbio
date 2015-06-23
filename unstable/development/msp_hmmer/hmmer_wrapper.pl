#!/usr/bin/env perl

use warnings;
use strict;
use File::Copy;

print "PATH: $ENV{'PATH'}";
my $cmd="@ARGV";
my ($tblout,$domtblout);
die("Missing --tblout\n") unless $cmd =~ / \-\-?tblout\s+(\S+)/;
$tblout = $1;
print "$tblout\n" ;
die("Missing --domtblout\n") unless $cmd =~ / \-\-?domtblout\s+(\S+)/;
$domtblout = $1;
print "$domtblout\n" ;
print "command= $cmd\n" ;
my $output = `$cmd 2>&1`;
die("HMMer failure: $output\n") unless $? == 0;
reformat_tblout($tblout);
reformat_domtblout($domtblout);
exit;

sub reformat_tblout
{
    my $infile = shift;
    my $outfile = "$infile.new";
    open(IN, "<$infile") or die($!);
    open(OUT, ">$outfile") or die($!);
    print OUT "#target name\taccession\tquery name\taccession\tE-value (full)\tscore (full)\tbias (full)\tE-value (best)\tscore (best)\tbias (best)\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription of target\n";
    while (<IN>)
    {
        next if /^#/;
        my @row0 = split(/\s+/);
        my @row = @row0[0..17];
        push @row, join(' ', @row0[18..$#row0]);
        print OUT join("\t", @row), "\n";
    }
    close(IN);
    close(OUT);
    move($outfile,$infile);
}

sub reformat_domtblout
{
    my $infile = shift;
    my $outfile = "$infile.new";
    open(IN, "<$infile") or die($!);
    open(OUT, ">$outfile") or die($!);
    print OUT "#target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value (full)\tscore (full)\tbias (full)\t# (this)\tof (this)\tc-Evalue (this)\ti-Evalue (this)\tscore (this)\tbias (this)\tfrom (hmm)\tto (hmm)\tfrom (ali)\tto (ali)\tfrom (env)\tto (env)\tacc\tdescription of target\n";
    while (<IN>)
    {
        next if /^#/;
        my @row0 = split(/\s+/);
        my @row = @row0[0..21];
        push @row, join(' ', @row0[22..$#row0]);
        print OUT join("\t", @row), "\n";
    }
    close(IN);
    close(OUT);
    move($outfile,$infile);
}
__END__
