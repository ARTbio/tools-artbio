#!/usr/bin/env python
# small RNA oriented bowtie wrapper
# version 1 19-5-2014
# Usage sRbowtie.py <1 input_fasta_file> <2 alignment method> <3 -v mismatches> <4 out_type> <5 buildIndexIfHistory> <6 fasta/bowtie index> <7 bowtie output> <8 ali_fasta> <9 unali_fasta> <10 --num-threads \${GALAXY_SLOTS:-4}>
# To Do:
# implement number of bowtie processes as a Galaxy env variable
# implement an arg parser
# Christophe Antoniewski <drosofff@gmail.com>

import sys, os, subprocess, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def bowtieCommandLiner (alignment_method, v_mis, out_type, aligned, unaligned, input, index, output, pslots="12"):
    if alignment_method=="RNA":
        x = "-v %s -M 1 --best --strata -p %s --norc --suppress 2,6,7,8" % (v_mis, pslots)
    elif alignment_method=="unique":
        x =  "-v %s -m 1 -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif  alignment_method=="multiple":
        x = "-v %s -M 1 --best --strata -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif alignment_method=="k_option":
        x = "-v %s -k 1 --best -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif alignment_method=="n_option":
        x = "-n %s -M 1 --best -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif alignment_method=="a_option":
        x = "-v %s -a --best -p %s --suppress 6,7,8" % (v_mis, pslots)
    if aligned == "None" and unaligned == "None": fasta_command = ""
    elif aligned != "None" and unaligned == "None": fasta_command= " --al %s" % aligned
    elif aligned == "None" and unaligned != "None": fasta_command = " --un %s" % unaligned
    else: fasta_command = " --al %s --un %s" % (aligned, unaligned)
    x = x + fasta_command
    if out_type == "tabular":
        return "bowtie %s %s -f %s > %s" % (x, index, input, output)
    elif out_type=="sam":
        return "bowtie %s -S %s -f %s > %s" % (x, index, input, output)
    elif out_type=="bam":
        return "bowtie %s -S %s -f %s |samtools view -bS - > %s" % (x, index, input, output)

def bowtie_squash(fasta):
  tmp_index_dir = tempfile.mkdtemp() # make temp directory for bowtie indexes
  ref_file = tempfile.NamedTemporaryFile( dir=tmp_index_dir )
  ref_file_name = ref_file.name
  ref_file.close() # by default, delete the temporary file, but ref_file.name is now stored in ref_file_name
  os.symlink( fasta, ref_file_name ) # symlink between the fasta source file and the deleted ref_file name
  cmd1 = 'bowtie-build -f %s %s' % (ref_file_name, ref_file_name ) # bowtie command line, which will work after changing dir (cwd=tmp_index_dir)
  try:
    FNULL = open(os.devnull, 'w')
    tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name # a path string for a temp file in tmp_index_dir. Just a string
    tmp_stderr = open( tmp, 'wb' ) # creates and open a file handler pointing to the temp file
    proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_index_dir, stderr=FNULL, stdout=FNULL ) # both stderr and stdout of bowtie-build are redirected in  dev/null
    returncode = proc.wait()
    tmp_stderr.close()
    FNULL.close()
    sys.stdout.write(cmd1 + "\n")
  except Exception, e:
    # clean up temp dir
    if os.path.exists( tmp_index_dir ):
      shutil.rmtree( tmp_index_dir )
      stop_err( 'Error indexing reference sequence\n' + str( e ) )
  # no Cleaning if no Exception, tmp_index_dir has to be cleaned after bowtie_alignment()
  index_full_path = os.path.join(tmp_index_dir, ref_file_name) # bowtie fashion path without extention
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
    first_command_line = " ".join(command_line.split()[:-3]) + " -o " + path_to_unsortedBam + " - "
    # example: bowtie -v 0 -M 1 --best --strata -p 12 --suppress 6,7,8 -S /home/galaxy/galaxy-dist/bowtie/Dmel/dmel-all-chromosome-r5.49 -f /home/galaxy/galaxy-dist/database/files/003/dataset_3460.dat |samtools view -bS -o /tmp/tmp_PgMT0/unsorted.bam -
    second_command_line = "samtools sort  %s %s" % (path_to_unsortedBam, path_to_sortedBam) # generates an "unsorted.bam.sorted.bam file", NOT an "unsorted.bam.sorted" file
    p = subprocess.Popen(args=first_command_line, cwd=tmp_index_dir, shell=True, stderr=tmp_stderr.fileno()) # fileno() method return the file descriptor number of tmp_stderr
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

def __main__():
  F = open (sys.argv[7], "w")
  if sys.argv[5] == "history":
    tmp_dir, index_path = bowtie_squash(sys.argv[6])
  else:
    tmp_dir, index_path = "dummy/dymmy", sys.argv[6]
  command_line = bowtieCommandLiner(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[8], sys.argv[9], sys.argv[1], index_path, sys.argv[7], sys.argv[10])
  bowtie_alignment(command_line, flyPreIndexed=tmp_dir)
  F.close()
if __name__=="__main__": __main__()
