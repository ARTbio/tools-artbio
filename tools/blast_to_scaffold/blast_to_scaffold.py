#!/usr/bin/env python
import sys
import argparse

def insert_newlines(string, every=60):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)
    
def getseq (fastadict, transcript, up, down, orientation="direct"):
    def reverse (seq):
        revdict = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
        revseq = [revdict[i] for i in seq[::-1]]
        return "".join(revseq)
    pickseq = fastadict[transcript][up-1:down]
    if orientation == "direct":
        return pickseq
    else:
        return reverse(pickseq)

def Parser():
    the_parser = argparse.ArgumentParser(
        description="Generate DNA scaffold from blastn or tblastx alignment of Contigs")
    the_parser.add_argument(
        '--sequences', action="store", type=str, help="input sequence file in fasta format")
    the_parser.add_argument(
        '--guideSequence', action="store", type=str, help="the reference sequence to guide the scaffold assembly in fasta format")
    the_parser.add_argument(
        '--blast-tab', dest="blast_tab", action="store", type=str, help="13-columns tabular blastn or tblastx output")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output file path, fasta format")
    args = the_parser.parse_args()
    return args
    
def blatnInfo (file):
    blastlist = []
    with open(file, "r") as f:
        for line in f:
            minilist = []
            fields = line.rstrip().split()
            minilist.append(fields[0])
            minilist.extend(fields[6:10])
            blastlist.append(minilist)
    blastlist.sort(key=lambda x: x[3], reverse=True)
    return blastlist
    
def myContigs (file):
    Contigs = {}
    with open(file, "r") as f:
        for line in f:
            if line[0] == ">":
                header = line[1:-1]
                Contigs[header] = ""
            else:
                Contigs[header] += line[:-1]
    return Contigs
    
def myGuide (file):
    Guide = {}
    coordinate = 0
    with open(file, "r") as f:
        for line in f:
            if line[0] == ">":
                continue
            else:
                for nucleotide in line[:-1]:
                    coordinate += 1
                    Guide[coordinate] = nucleotide.lower()
    return Guide

def updateGuide (blastlist, GuideDict, ContigsDict):
    '''
    the blastlist object is a list of list with
    element [0] : name of the blasted Contig
    element [1] : queryStart of the alignment to the reference
    element [2] = queryStop of the alignment to the reference        
    element [3] : subjectStart of the alignment to the reference
    element [4] = subjectStop of the alignment to the reference        
    '''
    for fields in blastlist:
        seqHeader = fields[0]
        queryStart = int(fields[1])
        queryStop = int(fields[2])
        subjectStart = int(fields[3])
        subjectStop = int(fields[4])
        if subjectStart > subjectStop:
            subjectStart, subjectStop = subjectStop, subjectStart
            orientation = "reverse"
        else:
            orientation = "direct"
        sequence = getseq (ContigsDict, seqHeader, queryStart, queryStop, orientation)
        for i in range(subjectStart, subjectStop+1):
            try:
                del GuideDict[i]
            except KeyError:
                continue
        for i, nucleotide in enumerate(sequence):
            GuideDict[i+subjectStart] = nucleotide
            
def finalAssembly (GuideDict, outputfile):
    finalSeqList = []
    for keys in sorted(GuideDict):
        finalSeqList.append(GuideDict[keys])
    finalSequence = insert_newlines("".join(finalSeqList) )
    Out = open (outputfile, "w")
    print >> Out, ">Scaffold"
    print >> Out, finalSequence
    Out.close()
    
def __main__():
    args = Parser()
    ContigsDict = myContigs (args.sequences)
    GuideDict = myGuide (args.guideSequence)
    blastlist = blatnInfo(args.blast_tab)
    updateGuide(blastlist, GuideDict, ContigsDict)
    finalAssembly(GuideDict, args.output)

if __name__ == "__main__":
    __main__()
