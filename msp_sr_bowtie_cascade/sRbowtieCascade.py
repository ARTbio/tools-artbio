#!/usr/bin/env python
# small RNA oriented bowtie wrapper in cascade for small RNA data set genome annotation
# version 0.9 13-6-2014
# Usage sRbowtie_cascade.py see Parser() for valid arguments
# Christophe Antoniewski <drosofff@gmail.com>

import sys, os, subprocess, tempfile, shutil, argparse
from collections import defaultdict

def Parser():
  the_parser = argparse.ArgumentParser()
  the_parser.add_argument('--output', action="store", type=str, help="output file")
  the_parser.add_argument('--num-threads', dest="num_threads", action="store", type=str, help="number of bowtie threads")
  the_parser.add_argument('--mismatch', action="store", type=str, help="number of mismatches allowed")
  the_parser.add_argument('--indexing-flags', dest="indexing_flags", nargs='+', help="whether the index should be generated or not by bowtie-buid")
  the_parser.add_argument('--index',nargs='+', help="paths to indexed or fasta references")
  the_parser.add_argument('--indexName',nargs='+', help="Names of the indexes")
  the_parser.add_argument('--input',nargs='+', help="paths to multiple input files")
  the_parser.add_argument('--label',nargs='+', help="labels of multiple input files")
  args = the_parser.parse_args()
  return args
 
def stop_err( msg ):
  sys.stderr.write( '%s\n' % msg )
  sys.exit()

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
  return index_full_path  
  
def make_working_dir():
  working_dir = tempfile.mkdtemp()
  return working_dir
  
def Clean_TempDir(directory):
  if os.path.exists( directory ):
    shutil.rmtree( directory )
  return

def bowtie_alignment(command_line="None", working_dir = ""):
  FNULL = open(os.devnull, 'w')
  p = subprocess.Popen(args=command_line, cwd=working_dir, shell=True, stderr=FNULL, stdout=FNULL)
  returncode = p.wait()
  sys.stdout.write("%s\n" % command_line)
  FNULL.close()
  #p = subprocess.Popen(["wc", "-l", "%s/al.fasta"%working_dir], cwd=working_dir, stdout=subprocess.PIPE)
  #aligned =  p.communicate()[0].split()[0]
  aligned = 0
  try: # hacked at gcc2014 in case of no alignment, no al.fasta file generated (?)
    F = open ("%s/al.fasta" % working_dir, "r")
    for line in F:
      aligned += 1
    F.close()
  except: pass
  sys.stdout.write("Aligned: %s\n" % aligned)
  return aligned/2

def CommandLiner (v_mis="1", pslots="12", index="dum/my", input="dum/my", working_dir=""):
  return "bowtie -v %s -k 1 --best -p %s --al %s/al.fasta --un %s/unal.fasta --suppress 1,2,3,4,5,6,7,8 %s -f %s" % (v_mis, pslots, working_dir, working_dir, index, input)

def __main__():
  args = Parser()
  ## first we make all indexes available. They can be already available or be squashed by bowtie-build
  ## we keep them in a list that alternates indexPath and "toClear" or "DoNotDelete"
  BowtieIndexList = []
  for indexing_flags, bowtiePath in zip (args.indexing_flags, args.index):
    if indexing_flags == "history":
      BowtieIndexList.append ( bowtie_squash (bowtiePath) )
      BowtieIndexList.append ( "toClear" )
    else:
      BowtieIndexList.append ( bowtiePath )
      BowtieIndexList.append ( "DoNotDelete") 
  ###### temporary Indexes are generated. They must be deleted at the end (after removing file name in the temp path) 
  ResultDict = defaultdict(list)
  for label, input in zip(args.label, args.input): ## the main cascade, iterating over samples and bowtie indexes
    workingDir = make_working_dir()
    cmd = CommandLiner (v_mis=args.mismatch, pslots=args.num_threads, index=BowtieIndexList[0], input=input, working_dir=workingDir)
    ResultDict[label].append( bowtie_alignment(command_line=cmd, working_dir = workingDir) ) # first step of the cascade
    if len(BowtieIndexList) > 2: # is there a second step to perform ?
      os.rename("%s/al.fasta"%workingDir, "%s/toAlign.fasta"%workingDir) ## end of first step. the aligned reads are the input of the next step
      cmd = CommandLiner (v_mis=args.mismatch, pslots=args.num_threads, index=BowtieIndexList[2], input="%s/toAlign.fasta"%workingDir, working_dir=workingDir)
      ResultDict[label].append( bowtie_alignment(command_line=cmd, working_dir = workingDir) )## second step of the cascade
    if len(BowtieIndexList) > 4:  ## remaining steps
      for BowtieIndexPath in BowtieIndexList[4::2]:
        try:
          os.unlink("%s/al.fasta" % workingDir) # hacked at gcc 2014, to remove previous al.fasta file that may interfere with counting if new al.fasta is empty
        except: pass
        os.rename("%s/unal.fasta"%workingDir, "%s/toAlign.fasta"%workingDir)
        cmd = CommandLiner (v_mis=args.mismatch, pslots=args.num_threads, index=BowtieIndexPath, input="%s/toAlign.fasta"%workingDir, working_dir=workingDir)
        ResultDict[label].append( bowtie_alignment(command_line=cmd, working_dir = workingDir) )
    Fun = open("%s/unal.fasta"%workingDir, "r") ## to finish, compute the number of unmatched reads
    n = 0
    for line in Fun:
      n += 1
    ResultDict[label].append(n/2)
    Fun.close()
    Clean_TempDir (workingDir) # clean the sample working directory
  ## cleaning
  for IndexPath, IndexFlag in zip(BowtieIndexList[::2], BowtieIndexList[1::2]):
    if IndexFlag == "toClear":
      Clean_TempDir ("/".join(IndexPath.split("/")[:-1]))
  ## end of cleaning
  
  
    
  F = open (args.output, "w")
  print >> F, "alignment reference\t%s" % "\t".join(args.label)
  for i, reference in enumerate(args.indexName):
    F.write ("%s" % reference)
    for sample in args.label:
      F.write ("\t%s" % "{:,}".format(ResultDict[sample][i]) )
    print >> F
  F.write ("Remaining Unmatched")
  for sample in args.label:
    F.write ("\t%s" % "{:,}".format(ResultDict[sample][-1]) ) 
  print >> F

  F.close()

if __name__=="__main__": __main__()
