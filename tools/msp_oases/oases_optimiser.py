#!/usr/bin/env python

"""
VelvetOptimiser Wrapper
refactored using the adaptation of
Konrad Paszkiewicz	University of Exeter, UK.

"""
import pkg_resources;
import logging, os, string, sys, tempfile, glob, shutil, types, urllib
import shlex, subprocess
from optparse import OptionParser, OptionGroup
from stat import *


def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    job_dir= os.getcwd()
    #tmp_work_dir = tempfile.mkdtemp(dir = job_dir) # make temp directory in the job_dir
    #Parse Command Line
    starthash = sys.argv[1]
    endhash = sys.argv[2]
    inputs = sys.argv[3]
    transcripts = sys.argv[4]
    transcripts_path = ''
    cmdline = "oases_optimiser.sh %s %s '%s' %s  2&1>/dev/null" % (starthash, endhash, inputs, job_dir) # 2&1>/dev/null
    print >> sys.stdout, cmdline # so will appear as blurb for file
    print >> sys.stdout, job_dir
    print >> sys.stdout, "PATH = %s" % (os.environ['PATH'])
    try:
        proc = subprocess.Popen( args=cmdline, shell=True, stderr=subprocess.PIPE ) #  cwd=job_dir
        returncode = proc.wait()
        # get stderr, allowing for case where it's very large
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += proc.stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        stop_err( 'Error running oases_optimiser.py' + str( e ) )
    out = open(transcripts,'w')
    transcript_path = os.path.join(job_dir, "MergedAssemblyFolder", 'transcripts.fa')
    print >> sys.stdout, transcript_path
    for line in open(transcript_path):
        out.write( "%s" % (line) )
    out.close()
 
if __name__ == "__main__": __main__()
