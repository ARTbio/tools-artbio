#!/usr/bin/env python

"""
VelvetOptimiser Wrapper
refactored using the adaptation of
Konrad Paszkiewicz	University of Exeter, UK.

"""
import os, sys
import subprocess


def stop_err(msg):
    sys.stderr.write("%s\n" % msg)
    sys.exit()


def oases_optimiser(starthash, endhash, input):
    '''
    Replaces call to oases_optimiser.sh. For all k-mers between
    starthash and endhash run velvet and oases.
    '''
    for i in xrange(starthash, endhash, 2):
        cmd1 = "velveth outputFolder_{0} {0} {1} && ".format(i, input)
        cmd2 = "velvetg outputFolder_{0} -read_trkg yes && ".format(i)
        cmd3 = "oases outputFolder_{0}".format(i)
        proc = subprocess.call(args=cmd1 + cmd2 + cmd3, shell=True, stdout=sys.stdout, stderr=sys.stdout)
        if not proc == 0:
            print("Oases failed at k-mer %s, skipping" % i)
            continue
    cmd4 = "velveth MergedAssemblyFolder 27 -long outputFolder_*/transcripts.fa && "
    cmd5 = "velvetg MergedAssemblyFolder -read_trkg yes -conserveLong yes && "
    cmd6 = "oases MergedAssemblyFolder -merge yes"
    proc = subprocess.call(args=cmd4 + cmd5 + cmd6, shell=True, stdout=sys.stdout, stderr=sys.stdout)
    if not proc == 0:
        raise Exception("Oases could not merge assembly")

def __main__():
    starthash = int(sys.argv[1])
    endhash = int(sys.argv[2])
    input = sys.argv[3]
    transcripts = sys.argv[4]
    try:
        oases_optimiser(starthash, endhash, input)
    except Exception, e:
        stop_err('Error running oases_optimiser.py\n' + str(e))
    with open(transcripts, 'w') as out:
        transcript_path = os.path.join("MergedAssemblyFolder", 'transcripts.fa')
        for line in open(transcript_path):
            out.write("%s" % (line))


if __name__ == "__main__": __main__()
