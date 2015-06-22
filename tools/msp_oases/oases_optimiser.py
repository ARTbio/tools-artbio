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

def oases_optimiser(starthash, endhash, input, job_dir):
    '''
    Replaces call to oases_optimiser.sh. For all k-mers between
    starthash and endhash run velvet and oases.
    '''
    for i in xrange(starthash, endhash, 2):
        cmd1="velveth {0}/outputFolder_{1} {1} {2} && ".format(job_dir, i, input)
        cmd2="velvetg {0}/outputFolder_{1} -read_trkg yes && ".format(job_dir, i)
        cmd3="oases {0}/outputFolder_{1}".format(job_dir, i)
        proc = subprocess.call( args=cmd1+cmd2+cmd3, shell=True, stdout=sys.stdout, stderr=sys.stderr )
    cmd4="velveth {0}/MergedAssemblyFolder 27 -long outputFolder_*/transcripts.fa && ".format(job_dir)
    cmd5="velvetg {0}/MergedAssemblyFolder -read_trkg yes -conserveLong yes && ".format(job_dir)
    cmd6="oases {0}/MergedAssemblyFolder -merge yes".format(job_dir)
    proc = subprocess.call( args=cmd4+cmd5+cmd6, shell=True, stdout=sys.stdout, stderr=sys.stderr )

def __main__():
    job_dir= os.getcwd()
    #Parse Command Line
    starthash = int(sys.argv[1])
    endhash = int(sys.argv[2])
    input = sys.argv[3]
    transcripts = sys.argv[4]
    transcripts_path = ''
    print >> sys.stdout, "PATH = %s" % (os.environ['PATH'])
    try:
        oases_optimiser(starthash, endhash, input, job_dir)
    except Exception, e:
        stop_err( 'Error running oases_optimiser.py' + str( e ) )
    out = open(transcripts,'w')
    transcript_path = os.path.join(job_dir, "MergedAssemblyFolder", 'transcripts.fa')
    print >> sys.stdout, transcript_path
    for line in open(transcript_path):
        out.write( "%s" % (line) )
    out.close()
 
if __name__ == "__main__": __main__()
