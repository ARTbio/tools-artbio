#!/usr/bin/env python
# mirdeep* python wrapper
# refactoring a python version of the MDS_wrapper.pl
# Usage MDS_wrapper.py <fasta input> <MDS genome> <gff3_output> <data.result> <data.cluster>

import sys, re, os, subprocess, shlex, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

MDS_path = "/home/galaxy/bin/MDS_command_line_v32/MDS_command_line"

input_full_path = sys.argv[1]
MDS_genome = sys.argv[2]
gff3_output = sys.argv[3]
dataresult = sys.argv[4]
datacluster = sys.argv[5]

tmp_MDS_work_dir = tempfile.mkdtemp(dir = MDS_path) # make temp directory for MDS analysis
os.symlink(input_full_path, tmp_MDS_work_dir+"/data.fa" ) # symlink between the fasta source file and the required "data.fa" input
os.symlink(MDS_path+"/MDS_command_line.jar", tmp_MDS_work_dir+"/MDS_command_line.jar" ) # symlink to jar source in working directory
os.symlink(MDS_path+"/genome", tmp_MDS_work_dir+"/genome")
os.symlink(MDS_path+"/targetScan", tmp_MDS_work_dir+"/targetScan")
os.symlink(MDS_path+"/targetScan_files", tmp_MDS_work_dir+"/targetScan_files")

# execute MirDeep*
command_line = "java -Xmx4g -jar " + MDS_path + "/MDS_command_line.jar -r 5 -g " + MDS_genome + " data.fa" # -Xmx12g
print command_line
#tmp = tempfile.NamedTemporaryFile( dir=tmp_MDS_work_dir ).name
#tmp_stderr = open( tmp, 'wb' )

try:
  os.chdir(tmp_MDS_work_dir)
  p = subprocess.Popen(args=command_line, cwd=tmp_MDS_work_dir, shell=True, stderr=sys.stderr)
  returncode = p.wait()
  shutil.copy2 ("data.result", dataresult)
  shutil.copy2 ("data.cluster", datacluster)
  dataFILE = open("data.result", "r")
  datafile = dataFILE.readlines()
  dataFILE.close()
  GFF3OUT = open(gff3_output, "w")
  print >> GFF3OUT,"##gff-version 3"
  print >> GFF3OUT, "##Seqid	Source	Type	Start	End	Score	Strand	Phase	Attributes"
  print >> GFF3OUT, "##"
  for line in datafile[1:]:
    fields = line.split("\t")
    Seqid, Source, Type, Start, End, Score, Strand, Phase = fields[2], "MirDeep*", "hairPin_loci", fields[4].split("-")[0], fields[4].split("-")[1], fields[1], fields[3], "."
    ID = "ID=%s;%s_reads;%s;%s;mature_seq:%s" % (fields[0],fields[5],fields[7],fields[8],fields[9])
    print >> GFF3OUT, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (Seqid, Source, Type, Start, End, Score, Strand, Phase, ID)
  GFF3OUT.close()
  if os.path.exists( tmp_MDS_work_dir ):
    shutil.rmtree( tmp_MDS_work_dir )
  else:
    print "Error in cleaning tmp working directory"

except Exception, e:
  # clean up temp dir
  if os.path.exists( tmp_MDS_work_dir ):
    shutil.rmtree( tmp_MDS_work_dir )
    stop_err( 'Error running MDS_command_line.jar\n' + str( e ) )


