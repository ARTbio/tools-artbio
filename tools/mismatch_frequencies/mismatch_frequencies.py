import pysam, re, string
import matplotlib.pyplot as plt
import pandas as pd
import json
from collections import defaultdict
from collections import OrderedDict
import argparse
import itertools

class MismatchFrequencies:
    '''Iterate over a SAM/BAM alignment file, collecting reads with mismatches. One
    class instance per alignment file. The result_dict attribute will contain a
    nested dictionary with name, readlength and mismatch count.'''
    def __init__(self, result_dict={}, alignment_file=None, name="name", minimal_readlength=21, 
                 maximal_readlength=21,
                 number_of_allowed_mismatches=1, 
                 ignore_5p_nucleotides=0, 
                 ignore_3p_nucleotides=0,
                 possible_mismatches = [
                        'AC', 'AG', 'AT',
                        'CA', 'CG', 'CT',
                        'GA', 'GC', 'GT',
                        'TA', 'TC', 'TG'
                ]):
    
        self.result_dict = result_dict
        self.name = name
        self.minimal_readlength = minimal_readlength
        self.maximal_readlength = maximal_readlength
        self.number_of_allowed_mismatches = number_of_allowed_mismatches
        self.ignore_5p_nucleotides = ignore_5p_nucleotides
        self.ignore_3p_nucleotides = ignore_3p_nucleotides
        self.possible_mismatches = possible_mismatches
        
        if alignment_file:
            self.pysam_alignment = pysam.Samfile(alignment_file)
            self.references = self.pysam_alignment.references #names of fasta reference sequences
            result_dict[name]=self.get_mismatches(
                self.pysam_alignment, 
                minimal_readlength, 
                maximal_readlength,
                possible_mismatches
            )
    
    def get_mismatches(self, pysam_alignment, minimal_readlength, 
                       maximal_readlength, possible_mismatches):
        mismatch_dict = defaultdict(int)
        rec_dd = lambda: defaultdict(rec_dd)
        len_dict = rec_dd()
        for alignedread in pysam_alignment:
            if self.read_is_valid(alignedread, minimal_readlength, maximal_readlength):
                chromosome = pysam_alignment.getrname(alignedread.rname)
                try:
                    len_dict[int(alignedread.rlen)][chromosome]['total valid reads'] += 1
                except TypeError:
                    len_dict[int(alignedread.rlen)][chromosome]['total valid reads'] = 1
                MD = alignedread.opt('MD')
                if self.read_has_mismatch(alignedread, self.number_of_allowed_mismatches):
                    (ref_base, mismatch_base)=self.read_to_reference_mismatch(MD, alignedread.seq, alignedread.is_reverse)
                    if ref_base == None:
                            continue
                    else:
                        for i, base in enumerate(ref_base):
                            if not ref_base[i]+mismatch_base[i] in possible_mismatches:
                                continue
                            try:
                                len_dict[int(alignedread.rlen)][chromosome][ref_base[i]+mismatch_base[i]] += 1
                            except TypeError:
                                len_dict[int(alignedread.rlen)][chromosome][ref_base[i]+mismatch_base[i]] = 1
        return len_dict
    
    def read_is_valid(self, read, min_readlength, max_readlength):
        '''Filter out reads that are unmatched, too short or
        too long or that contian insertions'''
        if read.is_unmapped:
            return False
        if read.rlen < min_readlength:
            return False
        if read.rlen > max_readlength:
            return False
        else:
            return True
    
    def read_has_mismatch(self, read, number_of_allowed_mismatches=1):
        '''keep only reads with one mismatch. Could be simplified'''
        NM=read.opt('NM')
        if NM <1: #filter out reads with no mismatch
            return False
        if NM >number_of_allowed_mismatches: #filter out reads with more than 1 mismtach
            return False
        else:
            return True
        
    def mismatch_in_allowed_region(self, readseq, mismatch_position):
        '''
        >>> M = MismatchFrequencies()
        >>> readseq = 'AAAAAA'
        >>> mismatch_position = 2
        >>> M.mismatch_in_allowed_region(readseq, mismatch_position)
        True
        >>> M = MismatchFrequencies(ignore_3p_nucleotides=2, ignore_5p_nucleotides=2)
        >>> readseq = 'AAAAAA'
        >>> mismatch_position = 1
        >>> M.mismatch_in_allowed_region(readseq, mismatch_position)
        False
        >>> readseq = 'AAAAAA'
        >>> mismatch_position = 4
        >>> M.mismatch_in_allowed_region(readseq, mismatch_position)
        False
        '''
        mismatch_position+=1 # To compensate for starting the count at 0
        five_p = self.ignore_5p_nucleotides
        three_p = self.ignore_3p_nucleotides
        if any([five_p > 0, three_p > 0]):
            if any([mismatch_position <= five_p, 
                    mismatch_position >= (len(readseq)+1-three_p)]): #Again compensate for starting the count at 0
                return False
            else:
                return True
        else:
            return True
            
    def read_to_reference_mismatch(self, MD, readseq, is_reverse):
        '''
        This is where the magic happens. The MD tag contains SNP and indel information,
        without looking to the genome sequence. This is a typical MD tag: 3C0G2A6.
        3 bases of the read align to the reference, followed by a mismatch, where the
        reference base is C, followed by 10 bases aligned to the reference. 
        suppose a reference 'CTTCGATAATCCTT'
                             |||  || ||||||
                 and a read 'CTTATATTATCCTT'. 
        This situation is represented by the above MD tag. 
        Given MD tag and read sequence this function returns the reference base C, G and A, 
        and the mismatched base A, T, T.
        >>> M = MismatchFrequencies()
        >>> MD='3C0G2A7'
        >>> seq='CTTATATTATCCTT'
        >>> result=M.read_to_reference_mismatch(MD, seq, is_reverse=False)
        >>> result[0]=="CGA"
        True
        >>> result[1]=="ATT"
        True
        >>> 
        '''
        search=re.finditer('[ATGC]',MD)
        if '^' in MD:
            print 'WARNING insertion detected, mismatch calling skipped for this read!!!'
            return (None, None)
        start_index=0 # refers to the leading integer of the MD string before an edited base
        current_position=0 # position of the mismatched nucleotide in the MD tag string
        mismatch_position=0 # position of edited base in current read 
        reference_base=""
        mismatched_base=""
        for result in search:
            current_position=result.start()
            mismatch_position=mismatch_position+1+int(MD[start_index:current_position]) #converts the leading characters before an edited base into integers
            start_index=result.end()
            reference_base+=MD[result.end()-1]
            mismatched_base+=readseq[mismatch_position-1]
        if is_reverse:
            reference_base=reverseComplement(reference_base)
            mismatched_base=reverseComplement(mismatched_base)
            mismatch_position=len(readseq)-mismatch_position-1
        if mismatched_base=='N':
            return (None, None)
        if self.mismatch_in_allowed_region(readseq, mismatch_position):
            return (reference_base, mismatched_base)
        else:
            return (None, None)

def reverseComplement(sequence):
    '''do a reverse complement of DNA base.
    >>> reverseComplement('ATGC')=='GCAT'
    True
    >>> 
    '''
    sequence=sequence.upper()
    complement = string.maketrans('ATCGN', 'TAGCN')
    return sequence.upper().translate(complement)[::-1]

def barplot(df, library, axes):
    df.plot(kind='bar', ax=axes, subplots=False,\
            stacked=False, legend='test',\
            title='Mismatch frequencies for {0}'.format(library))
    
def df_to_tab(df, output):
    df.to_csv(output, sep='\t')

def reduce_result(df, possible_mismatches):
    '''takes a pandas dataframe with full mismatch details and
    summarises the results for plotting.'''
    alignments = df['Alignment_file'].unique()
    readlengths = df['Readlength'].unique()
    combinations = itertools.product(*[alignments, readlengths]) #generate all possible combinations of readlength and alignment files
    reduced_dict = {}
    frames = []
    last_column = 3+len(possible_mismatches)
    for combination in combinations:
        library_subset = df[df['Alignment_file'] == combination[0]]
        library_readlength_subset = library_subset[library_subset['Readlength'] == combination[1]]
        sum_of_library_and_readlength = library_readlength_subset.iloc[:,3:last_column+1].sum()
        if not reduced_dict.has_key(combination[0]):
            reduced_dict[combination[0]] = {}
        reduced_dict[combination[0]][combination[1]] = sum_of_library_and_readlength.to_dict()
    return reduced_dict

def plot_result(reduced_dict, args):
    names=reduced_dict.keys()
    nrows=len(names)/2+1
    fig = plt.figure(figsize=(16,32))
    for i,library in enumerate (names):
        axes=fig.add_subplot(nrows,2,i+1)
        library_dict=reduced_dict[library]
        df=pd.DataFrame(library_dict)
        df.drop(['total aligned reads'], inplace=True)
        barplot(df, library, axes),
        axes.set_ylabel('Mismatch count / all valid reads * readlength')
    fig.savefig(args.output_pdf, format='pdf')    

def format_result_dict(result_dict, chromosomes, possible_mismatches):
    '''Turn nested dictionary into preformatted tab seperated lines'''
    header = "Reference sequence\tAlignment_file\tReadlength\t" + "\t".join(
        possible_mismatches) + "\ttotal aligned reads"
    libraries = result_dict.keys()
    readlengths = result_dict[libraries[0]].keys()
    result = []
    for chromosome in chromosomes:
        for library in libraries:
            for readlength in readlengths:
                line = []
                line.extend([chromosome, library, readlength])
                try:
                    line.extend([result_dict[library][readlength][chromosome].get(mismatch, 0) for mismatch in possible_mismatches])
                    line.extend([result_dict[library][readlength][chromosome].get(u'total valid reads', 0)])
                except KeyError:
                    line.extend([0 for mismatch in possible_mismatches])
                    line.extend([0])
                result.append(line)
    df = pd.DataFrame(result, columns=header.split('\t'))
    last_column=3+len(possible_mismatches)
    df['mismatches/per aligned nucleotides'] = df.iloc[:,3:last_column].sum(1)/(df.iloc[:,last_column]*df['Readlength'])
    return df
  
def setup_MismatchFrequencies(args):
    resultDict=OrderedDict()
    kw_list=[{'result_dict' : resultDict, 
             'alignment_file' :alignment_file, 
             'name' : name, 
             'minimal_readlength' : args.min, 
             'maximal_readlength' : args.max,
             'number_of_allowed_mismatches' : args.n_mm,
             'ignore_5p_nucleotides' : args.five_p, 
             'ignore_3p_nucleotides' : args.three_p,
             'possible_mismatches' : args.possible_mismatches }
             for alignment_file, name in zip(args.input, args.name)]
    return (kw_list, resultDict)

def nested_dict_to_df(dictionary):
    dictionary = {(outerKey, innerKey): values for outerKey, innerDict in dictionary.iteritems() for innerKey, values in innerDict.iteritems()}
    df=pd.DataFrame.from_dict(dictionary).transpose()
    df.index.names = ['Library', 'Readlength']
    return df

def run_MismatchFrequencies(args):
    kw_list, resultDict=setup_MismatchFrequencies(args)
    references = [MismatchFrequencies(**kw_dict).references for kw_dict in kw_list]
    return (resultDict, references[0])

def main():
    result_dict, references = run_MismatchFrequencies(args)
    df = format_result_dict(result_dict, references, args.possible_mismatches)
    reduced_dict = reduce_result(df, args.possible_mismatches)
    plot_result(reduced_dict, args)
    reduced_df = nested_dict_to_df(reduced_dict)
    df_to_tab(reduced_df, args.output_tab)
    if not args.expanded_output_tab == None:
        df_to_tab(df, args.expanded_output_tab)
    return reduced_dict

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Produce mismatch statistics for BAM/SAM alignment files.')
    parser.add_argument('--input', nargs='*', help='Input files in SAM/BAM format')
    parser.add_argument('--name', nargs='*', help='Name for input file to display in output file. Should have same length as the number of inputs')
    parser.add_argument('--output_pdf', help='Output filename for graph')
    parser.add_argument('--output_tab', help='Output filename for table')
    parser.add_argument('--expanded_output_tab', default=None, help='Output filename for table')
    parser.add_argument('--possible_mismatches', default=[
            'AC', 'AG', 'AT','CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG'
        ], nargs='+', help='specify mismatches that should be counted for the mismatch frequency. The format is Reference base -> observed base, eg AG for A to G mismatches.')
    parser.add_argument('--min', '--minimal_readlength', type=int, help='minimum readlength')
    parser.add_argument('--max', '--maximal_readlength', type=int, help='maximum readlength')
    parser.add_argument('--n_mm', '--number_allowed_mismatches', type=int, default=1, help='discard reads with more than n mismatches')
    parser.add_argument('--five_p', '--ignore_5p_nucleotides', type=int, default=0, help='when calculating nucleotide mismatch frequencies ignore the first N nucleotides of the read')
    parser.add_argument('--three_p', '--ignore_3p_nucleotides', type=int, default=1, help='when calculating nucleotide mismatch frequencies ignore the last N nucleotides of the read')
    #args = parser.parse_args(['--input', '3mismatches_ago2ip_s2.bam', '3mismatches_ago2ip_ovary.bam','--possible_mismatches','AC','AG', 'CG', 'TG', 'CT','--name', 'Siomi1', 'Siomi2' , '--five_p', '3','--three_p','3','--output_pdf', 'out.pdf', '--output_tab', 'out.tab', '--expanded_output_tab', 'expanded.tab', '--min', '20', '--max', '22'])
    args = parser.parse_args()
    reduced_dict = main()


