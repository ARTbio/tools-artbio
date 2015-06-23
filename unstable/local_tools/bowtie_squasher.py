#!/usr/bin/env python
# generic bowtie wrapper to make bowtie index
# version 2 (8-11-2012)
# Usage bowtie_squasher <fasta file> <input file> <output file> <-v option (0 to 3)>

import sys, re, os, subprocess, shlex, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def bowtie_squash(fasta, aligning_cmds, input_file, output_file, unaligned):
  stdout = ''
  # make temp directory for placement of indices and copy reference file there if necessary
  tmp_index_dir = tempfile.mkdtemp()
  ref_file = tempfile.NamedTemporaryFile( dir=tmp_index_dir )
  ref_file_name = ref_file.name
  ref_file.close()
  os.symlink( fasta, ref_file_name )
  cmd1 = 'bowtie-build -f %s %s' % (ref_file_name, ref_file_name )
  try:
    tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
    tmp_stderr = open( tmp, 'wb' )
    proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno() )
    returncode = proc.wait()
    tmp_stderr.close()
    # get stderr, allowing for case where it's very large
    tmp_stderr = open( tmp, 'rb' )
    stderr = ''
    buffsize = 1048576
    try:
      while True:
        stderr += tmp_stderr.read( buffsize )
        if not stderr or len( stderr ) % buffsize != 0:
          break
    except OverflowError:
      pass
    tmp_stderr.close()
    if returncode != 0:
      raise Exception, stderr
  except Exception, e:
    # clean up temp dir
    if os.path.exists( tmp_index_dir ):
      shutil.rmtree( tmp_index_dir )
      stop_err( 'Error indexing reference sequence\n' + str( e ) )
      stdout += 'File indexed. '
  # prepare actual mapping commands
  cmd2 = 'bowtie %s %s -f %s --un %s > %s' % ( aligning_cmds, ref_file_name, input_file, unaligned, output_file )
  # align
  tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
  tmp_stderr = open( tmp, 'wb' )
  proc = subprocess.Popen( args=cmd2, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno() )
  returncode = proc.wait()
  tmp_stderr.close()
  # get stderr, allowing for case where it's very large
  tmp_stderr = open( tmp, 'rb' )
  stderr = ''
  buffsize = 1048576
  try:
    # have to nest try-except in try-finally to handle 2.4
    try:
      while True:
        stderr += tmp_stderr.read( buffsize )
        if not stderr or len( stderr ) % buffsize != 0:
          break
    except OverflowError:
      pass
      tmp_stderr.close()
      if returncode != 0:
        raise Exception, stderr
      # check that there are results in the output file
      if os.path.getsize( options.output ) == 0:
        raise Exception, 'The output file is empty, there may be an error with your input file or settings.'
    except Exception, e:
      stop_err( 'Error aligning sequence. ' + str( e ) )
  finally:
    # clean up temp dir
    if os.path.exists( tmp_index_dir ):
      shutil.rmtree( tmp_index_dir )
  stdout += 'Sequence file aligned.\n'
  sys.stdout.write( stdout )

def __main__():
  command = "-v " + sys.argv[4] + " -M 1 --best --strata -p 12 --suppress 6,7,8"
  bowtie_squash(sys.argv[1], command, sys.argv[2], sys.argv[3], sys.argv[5])

if __name__=="__main__": __main__()

