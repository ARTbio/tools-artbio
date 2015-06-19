#!/usr/bin/python
#  blastn blastx parser revised debugged: 3-4-2015. Commit issue.
# drosofff@gmail.com

import sys
import argparse
from collections import defaultdict

def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--blast', action="store", type=str, help="Path to the blast output (tabular format, 12 column)")
    the_parser.add_argument('--sequences', action="store", type=str, help="Path to the fasta file with blasted sequences")
    the_parser.add_argument('--fastaOutput', action="store", type=str, help="fasta output file of blast hits")
    the_parser.add_argument('--tabularOutput', action="store", type=str, help="tabular output file of blast analysis")
    the_parser.add_argument('--flanking', action="store", type=int, help="number of flanking nucleotides added to the hit sequences") 
    the_parser.add_argument('--mode', action="store", choices=["verbose", "short"], type=str, help="reporting (verbose) or not reporting (short) oases contigs")
    args = the_parser.parse_args()
    if not all ( (args.sequences, args.blast, args.fastaOutput, args.tabularOutput) ):
        the_parser.error('argument(s) missing, call the -h option of the script')
    if not args.flanking:
        args.flanking = 0
    return args

def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    if len(lst) %2 == 0:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

def mean(lst):
    if len(lst) < 1:
        return 0
    return sum(lst) / float(len(lst))

def getfasta (fastafile):
    fastadic = {}
    for line in open (fastafile):
        if line[0] == ">":
            header = line[1:-1]
            fastadic[header] = ""
        else:
            fastadic[header] += line
    for header in fastadic:
        fastadic[header] = "".join(fastadic[header].split("\n"))
    return fastadic

def insert_newlines(string, every=60):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)
    
def getblast (blastfile):
    '''blastinfo [0]	Percentage of identical matches
       blastinfo [1]	Alignment length
       blastinfo [2]	Number of mismatches
       blastinfo [3]	Number of gap openings
       blastinfo [4]	Start of alignment in query
       blastinfo [5]	End of alignment in query
       blastinfo [6]	Start of alignment in subject (database hit)
       blastinfo [7]	End of alignment in subject (database hit)
       blastinfo [8]	Expectation value (E-value)
       blastinfo [9]	Bit score
       blastinfo [10]	Subject length (NEED TO BE SPECIFIED WHEN RUNNING BLAST) '''
    blastdic = defaultdict (dict) 
    for line in open (blastfile):
        fields = line[:-1].split("\t")
        transcript = fields[0]
        subject = fields[1]
        blastinfo = [float(fields[2]) ] # blastinfo[0]
        blastinfo = blastinfo + [int(i) for i in fields[3:10] ] # blastinfo[1:8] insets 1 to 7
        blastinfo.append(fields[10]) # blastinfo[8] E-value remains as a string type
        blastinfo.append(float(fields[11])) # blastinfo[9] Bit score
        blastinfo.append(int(fields[12])) # blastinfo[10] Subject length MUST BE RETRIEVED THROUGH A 13 COLUMN BLAST OUTPUT
        try:
            blastdic[subject][transcript].append(blastinfo)
        except:
            blastdic[subject][transcript] = [ blastinfo ]
    return blastdic

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

def subjectCoverage (fastadict, blastdict, subject, QueriesFlankingNucleotides=0):
    SubjectCoverageList = []
    HitDic = {}
    bitScores = []
    for transcript in blastdict[subject]:
        prefix = "%s--%s_" % (subject, transcript)
        hitNumber = 0
        for hit in blastdict[subject][transcript]:
            hitNumber += 1
            suffix = "hit%s_IdMatch=%s,AligLength=%s,E-val=%s" % (hitNumber, hit[0], hit[1], hit[8])
            HitDic[prefix+suffix] = GetHitSequence (fastadict, transcript, hit[4], hit[5], QueriesFlankingNucleotides) #query coverage by a hit is in hit[4:6]
            SubjectCoverageList += range (min([hit[6], hit[7]]), max([hit[6], hit[7]]) + 1) # subject coverage by a hit is in hit[6:8]
            bitScores.append(hit[9])
            subjectLength = hit [10] # always the same value for a given subject. Stupid but simple
    TotalSubjectCoverage = len ( set (SubjectCoverageList) )
    RelativeSubjectCoverage = TotalSubjectCoverage/float(subjectLength)
    return HitDic, subjectLength, TotalSubjectCoverage, RelativeSubjectCoverage, max(bitScores), mean(bitScores)
    
def GetHitSequence (fastadict, FastaHeader, leftCoordinate, rightCoordinate, FlankingValue):
    if rightCoordinate > leftCoordinate:
        polarity = "direct"
    else:
        polarity = "reverse"
        leftCoordinate, rightCoordinate = rightCoordinate, leftCoordinate
    if leftCoordinate - FlankingValue > 0:
        leftCoordinate -= FlankingValue
    else:
        leftCoordinate = 1
    return getseq (fastadict, FastaHeader, leftCoordinate, rightCoordinate, polarity)
    
def outputParsing (F, Fasta, results, Xblastdict, fastadict, mode="verbose"):
    F= open(F, "w")
    Fasta=open(Fasta, "w")
    if mode == "verbose":
        print >>F, "# SeqId\t%Identity\tAlignLength\tStartSubject\tEndSubject\t%QueryHitCov\tE-value\tBitScore\n"
        for subject in sorted (results, key=lambda x: results[x]["meanBitScores"], reverse=True):
            print >> F, "#\n# %s" % subject
            print >> F, "# Suject Length: %s" % (results[subject]["subjectLength"])
            print >> F, "# Total Subject Coverage: %s" % (results[subject]["TotalCoverage"])
            print >> F, "# Relative Subject Coverage: %s" % (results[subject]["RelativeSubjectCoverage"])
            print >> F, "# Maximum Bit Score: %s" % (results[subject]["maxBitScores"])
            print >> F, "# Mean Bit Score: %s" % (results[subject]["meanBitScores"])
            for header in results[subject]["HitDic"]:
                print >> Fasta, ">%s\n%s" % (header, insert_newlines(results[subject]["HitDic"][header]) )
            print >> Fasta, "" # final carriage return for the sequence
            for transcript in Xblastdict[subject]:
                transcriptSize = float(len(fastadict[transcript]))
                for hit in Xblastdict[subject][transcript]:
                    percentIdentity, alignLenght, subjectStart, subjectEnd, queryCov = hit[0], hit[1], hit[6], hit[7], "%.1f" % (abs(hit[5]-hit[4])/transcriptSize*100)
                    Eval, BitScore = hit[8], hit[9]
                    info = [transcript] + [percentIdentity, alignLenght, subjectStart, subjectEnd, queryCov, Eval, BitScore]
                    info = [str(i) for i in info]
                    info = "\t".join(info)
                    print >> F, info
    else:
        print >>F, "# subject\tsubject length\tTotal Subject Coverage\tRelative Subject Coverage\tMaximum Bit Score\tMean Bit Score"
        for subject in sorted (results, key=lambda x: results[x]["meanBitScores"], reverse=True):
            line = []
            line.append(subject)
            line.append(results[subject]["subjectLength"])
            line.append(results[subject]["TotalCoverage"])
            line.append(results[subject]["RelativeSubjectCoverage"])
            line.append(results[subject]["maxBitScores"])
            line.append(results[subject]["meanBitScores"])
            line = [str(i) for i in line]
            print >> F, "\t".join(line)
            for header in results[subject]["HitDic"]:
                print >> Fasta, ">%s\n%s" % (header, insert_newlines(results[subject]["HitDic"][header]) )
            print >> Fasta, "" # final carriage return for the sequence
    F.close()
    Fasta.close()
        
    

def __main__ ():
    args = Parser()
    fastadict = getfasta (args.sequences)
    Xblastdict = getblast (args.blast)
    results = defaultdict(dict)
    for subject in Xblastdict:
        results[subject]["HitDic"], results[subject]["subjectLength"], results[subject]["TotalCoverage"], results[subject]["RelativeSubjectCoverage"], results[subject]["maxBitScores"], results[subject]["meanBitScores"]  = subjectCoverage(fastadict, Xblastdict, subject, args.flanking)
    outputParsing (args.tabularOutput, args.fastaOutput, results, Xblastdict, fastadict, args.mode)
if __name__=="__main__": __main__()
