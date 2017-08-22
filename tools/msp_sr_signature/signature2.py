#!/usr/bin/python

import argparse

def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str, help="input alignment file")
    the_parser.add_argument(
        '--minquery', type=int, help="Minimum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxquery', type=int, help="Maximum readsize of query reads (nt) - must be an integer")
    the_parser.add_argument(
        '--mintarget', type=int, help="Minimum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--maxtarget', type=int, help="Maximum readsize of target reads (nt) - must be an integer")
    the_parser.add_argument(
        '--minscope', type=int, help="Minimum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--maxscope', type=int, help="Maximum overlap analyzed (nt) - must be an integer")
    the_parser.add_argument(
        '--outputOverlapDataframe', action="store", type=str, help="Overlap dataframe")
    the_parser.add_argument('--referenceGenome', action='store',
                            help="path to the bowtie-indexed or fasta reference")
    the_parser.add_argument('--extract_index', action='store_true',
                            help="specify if the reference is an indexed Bowtie reference")
    the_parser.add_argument('--graph', action='store', choices=[
                            "global", "lattice"], help="small RNA signature is computed either globally or by item (global-lattice)")
    the_parser.add_argument(
        '--rcode', type=str, help="R code to be passed to the python script")
    args = the_parser.parse_args()
    return args

args = Parser()

class Map:

    def __init__(self, bam_file, sample):
        self.sample_name = sample
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.chromosomes = dict(zip(self.bam_object.references,
                                self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object)

    def create_map(self, bam_object):
        '''
        Returns a map_dictionary {(chromosome,read_position,polarity):
                                                    [read_length, ...]}
        '''
        map_dictionary = defaultdict(list)
        # get empty value for start and end of each chromosome
        for chrom in self.chromosomes:
            map_dictionary[(chrom, 1, 'F')] = []
            map_dictionary[(chrom, self.chromosomes[chrom], 'F')] = []
        for chrom in self.chromosomes:
            for read in bam_object.fetch(chrom):
                positions = read.positions  # a list of covered positions
                for pos in positions:
                    if not map_dictionary[(chrom, pos+1, 'F')]:
                        map_dictionary[(chrom, pos+1, 'F')] = []
                if read.is_reverse:
                    map_dictionary[(chrom, positions[-1]+1,
                                    'R')].append(read.query_alignment_length)
                else:
                    map_dictionary[(chrom, positions[0]+1,
                                    'F')].append(read.query_alignment_length)
        return map_dictionary

  def signature (self, minquery, maxquery, mintarget, maxtarget, scope, zscore="no", upstream_coord=None, downstream_coord=None):
    ''' migration to memory saving by specifying possible subcoordinates
    see the readcount method for further discussion
    scope must be a python iterable; scope define the *relative* offset range to be computed'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    query_range = range (minquery, maxquery+1)
    target_range = range (mintarget, maxtarget+1)
    Query_table = {}
    Target_table = {}
    frequency_table = dict ([(i, 0) for i in scope])
    for offset in self.readDict:
      if (abs(offset) < upstream_coord or abs(offset) > downstream_coord): continue
      for size in self.readDict[offset]:
        if size in query_range:
          Query_table[offset] = Query_table.get(offset, 0) + 1
        if size in target_range:
          Target_table[offset] = Target_table.get(offset, 0) + 1
    for offset in Query_table:
      for i in scope:
        frequency_table[i] += min(Query_table[offset], Target_table.get(-offset -i +1, 0))
    if minquery==mintarget and maxquery==maxtarget: ## added to incorporate the division by 2 in the method (26/11/2013), see signature_options.py and lattice_signature.py
      frequency_table = dict([(i,frequency_table[i]/2) for i in frequency_table])
    if zscore == "yes":
      z_mean = mean(frequency_table.values() )
      z_std = std(frequency_table.values() )
      if z_std == 0:
        frequency_table = dict([(i,0) for i in frequency_table] )
      else:
        frequency_table = dict([(i, (frequency_table[i]- z_mean)/z_std) for i in frequency_table] )    
    return frequency_table

  def hannon_signature (self, minquery, maxquery, mintarget, maxtarget, scope, upstream_coord=None, downstream_coord=None):
    ''' migration to memory saving by specifying possible subcoordinates see the readcount method for further discussion
    note that scope must be an iterable (a list or a tuple), which specifies the relative offsets that will be computed'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    query_range = range (minquery, maxquery+1)
    target_range = range (mintarget, maxtarget+1)
    Query_table = {}
    Target_table = {}
    Total_Query_Numb = 0
    general_frequency_table = dict ([(i,0) for i in scope])
    ## filtering the appropriate reads for the study
    for offset in self.readDict:
      if (abs(offset) < upstream_coord or abs(offset) > downstream_coord): continue
      for size in self.readDict[offset]:
        if size in query_range:
          Query_table[offset] = Query_table.get(offset, 0) + 1
          Total_Query_Numb += 1
        if size in target_range:
          Target_table[offset] = Target_table.get(offset, 0) + 1
    for offset in Query_table:
      frequency_table = dict ([(i,0) for i in scope])
      number_of_targets = 0
      for i in scope:
        frequency_table[i] += Query_table[offset] *  Target_table.get(-offset -i +1, 0)
        number_of_targets += Target_table.get(-offset -i +1, 0)
      for i in scope:
        try:
          general_frequency_table[i] += (1. / number_of_targets / Total_Query_Numb) * frequency_table[i]
        except ZeroDivisionError :
          continue
    return general_frequency_table      


# Run the R script that is defined in the xml using the Rscript binary
# provided with R.
R_command = "Rscript " + args.rcode
process = subprocess.Popen(R_command.split())
