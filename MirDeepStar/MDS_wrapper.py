#!/usr/bin/env python
# mirdeep* python wrapper
# refactoring a python version of the MDS_wrapper.pl
# Usage MDS_wrapper.py <fasta input> <MDS genome> <gff3_output>

import sys, re, os, subprocess, shlex, tempfile, shutil

MDS_path = "/home/galaxy/bin/MDS_command_line_v32/MDS_command_line"

input_full_path = sys.argv[1]
MDS_genome = sys.argv[2]
gff3_output = sys.argv[3]

# delete if previous data.fa and creates a symbolic link of input data
if os.path.isfile(path_to_sortedBam + ".bam"):
  shutil.copy2(path_to_sortedBam + ".bam", target_file)


if (-e "data.fa") {
  unlink("data.fa");
  };
symlink($input_full_path, $MDS_path."/data.fa");


# go to the MDS base directory
chdir($MDS_path) or die "cannot change: $!\n";

# execute MirDeep*
print "java -jar -Xmx12g MDS_command_line.jar -r 50 -g $MDS_genome data.fa\n" ;
`java -jar -Xmx12g MDS_command_line.jar -r 50 -g $MDS_genome data.fa`;

# parse data.result with MDSstar2gff3.py to generate gff3 results

`MDSstar2gff3.py data.result > $gff3_output `;

unlink ("data.fa");


def bowtie_squash(fasta):
  # make temp directory for placement of indices and copy reference file there if necessary
  tmp_index_dir = tempfile.mkdtemp()
  ref_file = tempfile.NamedTemporaryFile( dir=tmp_index_dir )
  ref_file_name = ref_file.name
  ref_file.close() # by default, this action delete the temporary file !
  os.symlink( fasta, ref_file_name ) # now there is a symlink between the fasta source file and the deleted ref_file name
  cmd1 = 'bowtie-build -f %s %s' % (ref_file_name, ref_file_name ) # this will work with the bowtie command line but we have to change the working dir !
  try:
    FNULL = open(os.devnull, 'w')
    tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name # a full path temp file to tmp_index_dir is created. just a string
    tmp_stderr = open( tmp, 'wb' ) # creates and open the file handler
    proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_index_dir, stderr=FNULL, stdout=FNULL ) # fileno() method return the file descriptor number of tmp_stderr / stderr=tmp_stderr.fileno()
    # here I played a while a finish by redirecting everythin in dev/null. Clean later tmp_stderr calls
    returncode = proc.wait()
    tmp_stderr.close()
    FNULL.close()
    sys.stdout.write(cmd1 + "\n")
  except Exception, e:
    # clean up temp dir
    if os.path.exists( tmp_index_dir ):
      shutil.rmtree( tmp_index_dir )
      stop_err( 'Error indexing reference sequence\n' + str( e ) )
  # no Cleaning if no Exception, to be cleaned later after bowtie alignment
  index_full_path = os.path.join(tmp_index_dir, ref_file_name) # bowtie fashion path (without extention) ...
  return tmp_index_dir, index_full_path

def bowtie_alignment(command_line, flyPreIndexed=''):
  # make temp directory just for stderr
  tmp_index_dir = tempfile.mkdtemp()
  tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
  tmp_stderr = open( tmp, 'wb' )
  # conditional statement for sorted bam generation viewable in Trackster
  if "samtools" in command_line:
    target_file = command_line.split()[-1] # recover the final output file name
    path_to_unsortedBam = os.path.join(tmp_index_dir, "unsorted.bam")
    path_to_sortedBam = os.path.join(tmp_index_dir, "unsorted.bam.sorted")
    first_command_line = " ".join(command_line.split()[:-3]) + " -o " + path_to_unsortedBam + " - " # example: bowtie -v 0 -M 1 --best --strata -p 12 --suppress 6,7,8 -S /home/galaxy/galaxy-dist/bowtie/Dmel/dmel-all-chromosome-r5.49 -f /home/galaxy/galaxy-di$
    second_command_line = "samtools sort  %s %s" % (path_to_unsortedBam, path_to_sortedBam) # Be carreful : this indeed will generate an unsorted.bam.sorted.bam file, NOT a unsorted.bam.sorted file
    p = subprocess.Popen(args=first_command_line, cwd=tmp_index_dir, shell=True, stderr=tmp_stderr.fileno())
    returncode = p.wait()
    sys.stdout.write("%s\n" % first_command_line + str(returncode))
    p = subprocess.Popen(args=second_command_line, cwd=tmp_index_dir, shell=True, stderr=tmp_stderr.fileno())
    returncode = p.wait()
    sys.stdout.write("\n%s\n" % second_command_line + str(returncode))
    if os.path.isfile(path_to_sortedBam + ".bam"):
      shutil.copy2(path_to_sortedBam + ".bam", target_file)
  else:
    p = subprocess.Popen(args=command_line, shell=True, stderr=tmp_stderr.fileno())
    returncode = p.wait()
    sys.stdout.write(command_line + "\n")
  tmp_stderr.close()
  ## cleaning if the index was created in the fly
  if os.path.exists( flyPreIndexed ):
    shutil.rmtree( flyPreIndexed )
  # cleaning tmp files and directories
  if os.path.exists( tmp_index_dir ):
    shutil.rmtree( tmp_index_dir )
  return

