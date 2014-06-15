#!/usr/bin/env python
# small RNA oriented bowtie wrapper in cascade for small RNA data set genome annotation
# version 0.9 13-6-2014
# Usage sRbowtie_cascade.py <1 analysis_output> <2 --num-threads \${GALAXY_SLOTS:-4}> <3 -v mismatches>
#                           <4 fastaInputs,,,> <5 buildIndexIfHistory,fasta/bowtieIndex,ad lib...>
# Christophe Antoniewski <drosofff@gmail.com>

import sys, os, subprocess, tempfile, shutil, argparse, collections

def Parser():
  the_parser = argparse.ArgumentParser()
  the_parser.add_argument('--output', action="store", type=str, help="output file")
  the_parser.add_argument('--num-threads', dest="--num_threads", action="store", type=str, help="number of bowtie threads")
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

def CommandLiner (v_mis="1", pslots="12", index="dum/my", input="dum/my"):
  return "bowtie -v %s -k 1 --best -p %s --al al.fasta --un unal.fasta --suppress 1,2,3,4,5,6,7,8 %s -f %s" % (v_mis, pslots, index, input)

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
  sys.stdout.write("%s\n" % command_line + str(returncode))
  FNULL.close()
  p = subprocess.Popen(["wc", "-l", "al.fasta"], cwd=working_dir, stdout=subprocess.PIPE)
  aligned =  p.communicate()[0].split()[0]
  return int(aligned)/2

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
      BowtieIndexList.append ( "DoNotDelete") )
  ###### temporary Indexes are generated. They must be deleted at the end (after removing file name in the temp path) 
  workingDir = make_working_dir()
  ResultDict = defaultdict(list)
  F = open (args.output, "w")
  for label, input in zip(args.label, args.input): ## the main cascade, iterating over samples and bowtie indexes
    cmd = CommandLiner (v_mis=args.mismatch, pslots=args.num_threads, index=BowtieIndexList[0], input=input)
    ResultDict[label].append( bowtie_alignment(command_line="cmd", working_dir = "") )
    cmd = CommandLiner (v_mis=args.mismatch, pslots=args.num_threads, index=BowtieIndexList[2], input=al.fasta)
    ResultDict[label].append( bowtie_alignment(command_line="cmd", working_dir = "") )
    if len(BowtieIndexList) > 4:
      for BowtieIndexPath in BowtieIndexList[4::2]:
        cmd = CommandLiner (v_mis=args.mismatch, pslots=args.num_threads, index=BowtieIndexPath, input="unal.fasta")
        ResultDict[label].append( bowtie_alignment(command_line="None", working_dir = "") )
  ## cleaning
  for IndexPath in BowtieIndexList[::2]:
    Clean_TempDir ("/".join(IndexPath.split("/")[:-1]))
  Clean_TempDir (workingDir)
  ## end of cleaning
  
  
    
  F = open (args.output, "w")

if __name__=="__main__": __main__()