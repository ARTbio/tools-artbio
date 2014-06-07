#!/usr/bin/perl
# my poor mirdeep* perl wrapper
# Usage MDS_wrapper.pl <fasta input> <MDS genome> <gff3_output> <data.result> <data.cluster>

use Cwd;
use File::Path 'rmtree';

my $original = cwd;
my $timestamp = time();
my $tempdir = $timestamp;

$MDS_path = "/home/galaxy/bin/MDS_command_line_v32/MDS_command_line";

$input_full_path = shift @ARGV ;
$MDS_genome = shift @ARGV ;
$gff3_output = shift @ARGV ;
$dataresult = shift @ARGV ;
$datacluster = shift @ARGV ;

# create and go to the MDS temp directory
mkdir $MDS_path.'/'.$tempdir or die "cannot create temp dir\n";
chdir($MDS_path.'/'.$tempdir) or die "cannot change: $!\n";

# creates a symbolic link of input data and MDS environment
symlink($input_full_path, "data.fa");
symlink("../MDS_command_line.jar", "MDS_command_line.jar");
symlink("../genome", "genome");
symlink("../targetScan", "targetScan");
symlink("../targetScan_files", "targetScan_files");

# execute MirDeep*
print "java -jar -Xmx12g MDS_command_line.jar -r 5 -g $MDS_genome data.fa\n" ;
print "in temp directory $MDS_path /$tempdir";
`java -jar -Xmx12g MDS_command_line.jar -r 5 -g $MDS_genome data.fa`;

# parse data.result with MDSstar2gff3.py to generate gff3 results

`MDSstar2gff3.py data.result > $gff3_output `;
`cp data.result $dataresult ` ;
`cp data.cluster $datacluster ` ;

chdir ("..");
rmtree ($tempdir) or die "cannot unlink the temp directory\n";

