#!/usr/bin/perl

use File::Copy qw/ copy /;

$fastafile = $ARGV[0];
$file_path = $ARGV[1];
$output_file = $ARGV[2];
$fastafile_text = $ARGV[3];

$output_basename = `basename $output_file`;
chomp $output_basename;
$filepath_basename = `basename $file_path`;
chomp $filepath_basename;

$output_dir = $output_file;
$output_dir =~ s/$output_basename/$filepath_basename/;

system ("mkdir -p $file_path");
copy ($fastafile, $output_file);
copy ($output_file, $file_path);

system ("bowtie-build $fastafile $file_path/$output_basename");

open (OUTPUT,">$output_file");
print OUTPUT "<h1>Bowtie index on $fastafile_text</h1>\n";
$dirout = `ls $file_path`;

foreach $file (split (/\n/, $dirout)) {
	print OUTPUT "<a href='$file'>$file</a><br/>\n";
}
close (OUTPUT);
