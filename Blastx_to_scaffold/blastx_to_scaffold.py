#!/usr/bin/python
import sys
import argparse

def insert_newlines(string, every=60):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)
    
def Parser():
    the_parser = argparse.ArgumentParser(
        description="Generate DNA scaffold from blastx alignment of Contigs")
    the_parser.add_argument(
        '--sequences', action="store", type=str, help="input sequence file in fasta format")
    the_parser.add_argument(
        '--blastx-tab', dest="blastx_tab", action="store", type=str, help="13-columns tabular blastx output")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output file path, fasta format")
    args = the_parser.parse_args()
    return args

def __main__():
    args = Parser()
    protLenght = int (open (args.blastx_tab, "r").readline().split("\t")[12])
    BlastxOutput = open (args.blastx_tab, "r")
    Contigs = open (args.sequences, "r")
    ContigsDict = {}
    protScaffold = {}
    
    for line in Contigs:
        if line[0] == ">":
            header = line[1:-1]
            ContigsDict[header] = ""
        else:
            ContigsDict[header] += line[:-1]
            
    protScaffold = dict ( [(i,"NNN") for i in range (1, protLenght+1)] )
        
    for line in BlastxOutput:
        fields = line[:-1].split("\t")
        queryStart = int(fields[6])
        queryStop = int(fields[7])
        subjectStart = int(fields[8])
        subjectStop = int(fields[9])
        seqHeader = fields[0]
        sequence = ContigsDict[seqHeader]
        for i in range (subjectStart, subjectStop):
            del protScaffold[i]
        protScaffold[subjectStop] = ContigsDict[seqHeader][queryStart -1: queryStop]
        
    finalSeqList = []
    for i in sorted (protScaffold):
        finalSeqList.append(protScaffold[i])
    finalSequence = insert_newlines("".join(finalSeqList))
     
    Out = open (args.output, "w")
    print >> Out, ">Scaffold"
    print >> Out, finalSequence
            
    BlastxOutput.close()
    Contigs.close()
    Out.close()
        
if __name__ == "__main__":
    __main__()
