#!/usr/bin/python
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
        description="Generate DNA scaffold from blastn alignment of Contigs")
    the_parser.add_argument(
        '--sequences', action="store", type=str, help="input sequence file in fasta format")
    the_parser.add_argument(
        '--guideSequence', action="store", type=str, help="the reference sequence to guide the scaffold assembly in fasta format")
    the_parser.add_argument(
        '--blastn-tab', dest="blastn_tab", action="store", type=str, help="13-columns tabular blastn output")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output file path, fasta format")
    args = the_parser.parse_args()
    return args

def __main__():
    args = Parser()
    dnaLenght = int (open (args.blastn_tab, "r").readline().split("\t")[12])
    BlastnOutput = open (args.blastn_tab, "r")
    Contigs = open (args.sequences, "r")
    ContigsDict = {}
    GuideDict = {}
    protScaffold = {}
    
    for line in Contigs:
        if line[0] == ">":
            header = line[1:-1]
            ContigsDict[header] = ""
        else:
            ContigsDict[header] += line[:-1]
            
    dnaScaffold = {}
    
    Guide = open(args.guideSequence, "r")
    coordinate = 0
    for line in Guide:
        if line[0] == ">":
            continue
        else:
            for nucleotide in line[:-1]:
                coordinate += 1
                GuideDict[coordinate] = nucleotide.lower()
    Guide.close()
    
        
    for line in BlastnOutput:
        fields = line[:-1].split("\t")
        queryStart = int(fields[6])
        queryStop = int(fields[7])
        subjectStart = int(fields[8])
        subjectStop = int(fields[9])
        if subjectStart > subjectStop:
            subjectStart, subjectStop = subjectStop, subjectStart
            orientation = "reverse"
        else:
            orientation = "direct"
        seqHeader = fields[0]
        sequence = getseq (ContigsDict, seqHeader, queryStart, queryStop, orientation)
        for i, nucleotide in enumerate(sequence):
            dnaScaffold[i+subjectStart] = nucleotide
        
    finalSeqList = []
    for i in range(1, dnaLenght+1) :
        try:
            finalSeqList.append(dnaScaffold[i])
        except KeyError:
            finalSeqList.append(GuideDict[i])
    finalSequence = insert_newlines("".join(finalSeqList))
     
    Out = open (args.output, "w")
    print >> Out, ">Scaffold"
    print >> Out, finalSequence
            
    BlastnOutput.close()
    Contigs.close()
    Out.close()
        
if __name__ == "__main__":
    __main__()
