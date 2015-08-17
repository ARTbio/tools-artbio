#!/usr/bin/perl

use File::Copy qw/ copy /;
use Cwd;

$mirna_file = $ARGV[0];
$trans_file = $ARGV[1];
$targets = $ARGV[2];
$file_path = $ARGV[3];

# do all the dirty work in a temp directory
$cwd = cwd();
system ("mkdir -p $file_path");
chdir($file_path);

$ret_ts = `TargetSpy -microRNAs $mirna_file -transcripts $trans_file -result targets`;
if ($ret_ts ne "") {die "TargetSpy error"}

system ("gunzip targets.gz");
copy ("targets", $targets);

chdir($cwd);
