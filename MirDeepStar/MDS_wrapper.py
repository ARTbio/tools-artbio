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
command_line = "java -jar -Xmx6g MDS_command_line.jar -r 5 -g " + "MDS_genome " + "data.fa"
tmp = tempfile.NamedTemporaryFile( dir=tmp_MDS_work_dir ).name
tmp_stderr = open( tmp, 'wb' )
try:
  p = subprocess.Popen(args=command_line, cwd=tmp_MDS_work_dir, stderr=tmp_stderr.fileno())
  returncode = p.wait()
except Exception, e:
  pass
  # clean up temp dir
  #if os.path.exists( tmp_MDS_work_dir ):
    #shutil.rmtree( tmp_MDS_work_dir )
    #stop_err( 'Error running MDS_command_line.jar\n' + str( e ) )

