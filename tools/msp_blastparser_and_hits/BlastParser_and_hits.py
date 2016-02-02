#!/usr/bin/python
#  blastn tblastn blastx parser revised 14-1-2016.
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
    the_parser.add_argument('--filter_relativeCov', action="store", type=float, default=0, help="filter out relative coverages below the specified ratio (float number)")
    the_parser.add_argument('--filter_maxScore', action="store", type=float, default=0, help="filter out best BitScores below the specified float number")
    the_parser.add_argument('--filter_meanScore', action="store", type=float, default=0, help="filter out mean BitScores below the specified float number")
    the_parser.add_argument('--filter_term_in', action="store", type=str, default="", help="select the specified term in the subject list")
    the_parser.add_argument('--filter_term_out', action="store", type=str, default="", help="exclude the specified term from the subject list")
    the_parser.add_argument('--al_sequences', action="store", type=str, help="sequences that have been blast aligned")
    the_parser.add_argument('--un_sequences', action="store", type=str, help="sequences that have not been blast aligned")
    the_parser.add_argument('--dataset_name', action="store", type=str, default="", help="the name of the dataset that has been parsed, to be reported in the output")
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
    
def outputParsing (dataset_name, F, Fasta, results, Xblastdict, fastadict, filter_relativeCov=0, filter_maxScore=0, filter_meanScore=0, filter_term_in="", filter_term_out="", mode="verbose"):
    def filter_results (results, filter_relativeCov=0, filter_maxScore=0, filter_meanScore=0, filter_term_in="", filter_term_out=""):
        for subject in results.keys():
            if results[subject]["RelativeSubjectCoverage"]<filter_relativeCov:
                del results[subject]
                continue
            if results[subject]["maxBitScores"]<filter_maxScore:
                del results[subject]
                continue
            if results[subject]["meanBitScores"]<filter_meanScore:
                del results[subject]
                continue
            if filter_term_in in subject:
                pass
            else:
                del results[subject]
                continue
            if filter_term_out and filter_term_out in subject:
                del results[subject]
                continue
        return results
        
    F= open(F, "w")
    Fasta=open(Fasta, "w")
    blasted_transcripts = []
    filter_results (results, filter_relativeCov, filter_maxScore, filter_meanScore, filter_term_in, filter_term_out)
    for subject in results:
        for transcript in Xblastdict[subject]:
            blasted_transcripts.append(transcript)
    blasted_transcripts = list( set( blasted_transcripts))
    if mode == "verbose":
        print >>F, "--- %s ---" % (dataset_name)
        print >>F, "# SeqId\t%Identity\tAlignLength\tStartSubject\tEndSubject\t%QueryHitCov\tE-value\tBitScore"
        for subject in sorted (results, key=lambda x: results[x]["meanBitScores"], reverse=True):
            print >> F, " \n# %s" % subject
            print >> F, "# Suject Length: %s" % (results[subject]["subjectLength"])
            print >> F, "# Total Subject Coverage: %s" % (results[subject]["TotalCoverage"])
            print >> F, "# Relative Subject Coverage: %s" % (results[subject]["RelativeSubjectCoverage"])
            print >> F, "# Best Bit Score: %s" % (results[subject]["maxBitScores"])
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
        print >>F, "--- %s ---" % (dataset_name)
        print >>F, "# subject\tsubject length\tTotal Subject Coverage\tRelative Subject Coverage\tBest Bit Score\tMean Bit Score"
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
    return blasted_transcripts
        
def dispatch_sequences (fastadict, blasted_transcripts, matched_sequences, unmatched_sequences):
    '''to output the sequences that matched and did not matched in the blast'''
    F_matched = open (matched_sequences, "w")
    F_unmatched = open (unmatched_sequences, "w")
    for transcript in fastadict:
        if transcript in blasted_transcripts: # le list of blasted_transcripts is generated by the outputParsing function
            print >> F_matched, ">%s\n%s" % (transcript, insert_newlines(fastadict[transcript]) )
        else:
            print >> F_unmatched, ">%s\n%s" % (transcript, insert_newlines(fastadict[transcript]) )
    F_matched.close()
    F_unmatched.close()
    return

def __main__ ():
    args = Parser()
    fastadict = getfasta (args.sequences)
    Xblastdict = getblast (args.blast)
    results = defaultdict(dict)
    for subject in Xblastdict:
        results[subject]["HitDic"], results[subject]["subjectLength"], results[subject]["TotalCoverage"], results[subject]["RelativeSubjectCoverage"], results[subject]["maxBitScores"], results[subject]["meanBitScores"]  = subjectCoverage(fastadict, Xblastdict, subject, args.flanking)
    blasted_transcripts = outputParsing (args.dataset_name, args.tabularOutput, args.fastaOutput, results, Xblastdict, fastadict,
                                        filter_relativeCov=args.filter_relativeCov, filter_maxScore=args.filter_maxScore,
                                        filter_meanScore=args.filter_meanScore, filter_term_in=args.filter_term_in,
                                        filter_term_out=args.filter_term_out, mode=args.mode)
    dispatch_sequences (fastadict, blasted_transcripts, args.al_sequences, args.un_sequences)

if __name__=="__main__": __main__()