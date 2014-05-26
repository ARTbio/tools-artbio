#!/usr/bin/env python
# Generic bowtie wrapper with -n option and bowtie output with column 2 (strand polarity)
# version 1
# Usage bowtie_rnaseq_wrapper.py <trimmed .rawfasta file> <index> <output>

import sys, re, os, subprocess, shlex, tempfile, shutil


def bowtie_alignment(rawfasta, index, output, bowtie_opt="-n 2 -p 12"):
  # make temp directory for placement of indices and copy reference file there if necessary
  tmp_index_dir = tempfile.mkdtemp()
  
  tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
  tmp_stderr = open( tmp, 'wb' )
  command_line = "bowtie %s %s -f %s > %s" % (bowtie_opt, index, rawfasta, output)
  print "#", command_line
  p = subprocess.Popen(args=command_line, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno()) #, stdout=bowtie_output
  returncode = p.wait()
  tmp_stderr.close()
  if os.path.exists( tmp_index_dir ):
    shutil.rmtree( tmp_index_dir )
  return
  
def __main__():
  bowtie_alignment(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__=="__main__": __main__()
