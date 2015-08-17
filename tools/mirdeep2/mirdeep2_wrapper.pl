#!/usr/bin/perl

use Cwd;

$collapsed_reads = $ARGV[0];
$bowtie_index_name = $ARGV[1];
$bowtie_files_path = $ARGV[2];
$arf_file = $ARGV[3];
$mirna_ref = $ARGV[4];
$mirna_other = $ARGV[5];
$precursors = $ARGV[6];
$file_path = $ARGV[7];
$input_name = $ARGV[8];

$csv_output = $ARGV[9];
$html_output = $ARGV[10];
$survey_output = $ARGV[11];
$mrd_output = $ARGV[12];

# the rest are options
$options = join (" ", @ARGV[13..$#ARGV]);

#point to correct bowtie index path
$basename = `basename $bowtie_index_name`;
chomp $basename;

# create the path used to house the pdfs
system ("mkdir -p $file_path");

# do all the dirty work in a temp directory
$cwd = cwd();
chdir($file_path);

if ($mirna_ref eq "None") {$mirna_ref="none";}
if ($mirna_other eq "None") {$mirna_other="none";}
if ($precursors eq "None") {$precursors="none";}

$ret_mirdeep2 = `miRDeep2.pl $collapsed_reads $bowtie_files_path/$basename $arf_file $mirna_ref $mirna_other $precursors -v $options 2> /dev/null`;
system ("mv $file_path/pdfs*/*.pdf $file_path");

# if ($csv_file eq "") {die "Error: Cannot find csv result file in dir $file_path\n";}
# if ($html_file eq "") {die "Error: Cannot find html result file in dir $file_path\n";}
# if ($survey_file eq "") {die "Error: Cannot find survey result file in dir $file_path\n";}
# if ($mrd_file eq "") {die "Error: Cannot find hairpin result file in dir $file_path\n";}

# replacing mirdeep created links with relative links in galaxy
system ("sed -r -i 's/file:\\/\\/.+\\/(.+\\.pdf)/\\1/g' $file_path/result*.html");

# hack to deal with galaxy html tooltip rendering issues
$html = `cat $file_path/result*.html`;
$html =~ s/(class=\"tooltip\d*\")>(.+)<span>(.+)<\/span>/\1 title=\"\3\">\2/g;
$html_file = `ls $file_path/result*.html`;
chomp $html_file;
open (HTML, ">$html_file");
print HTML $html;
close (HTML);

system ("cp -f $file_path/result*.csv $csv_output");
system ("cp -f $file_path/result*.html $html_output");
system ("cp -f $file_path/mirdeep_runs/run*/survey.csv $survey_output");
system ("cp -f $file_path/mirdeep_runs/run*/output.mrd $mrd_output");

chdir($cwd);
