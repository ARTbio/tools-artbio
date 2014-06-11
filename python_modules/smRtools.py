#!/usr/bin/python
# version 1 7-5-2012 unification of the SmRNAwindow class

import sys, subprocess
from collections import defaultdict, OrderedDict
from numpy import mean, median, std
from scipy import stats

def get_fasta (index="/home/galaxy/galaxy-dist/bowtie/5.37_Dmel/5.37_Dmel"):
  '''This function will return a dictionary containing fasta identifiers as keys and the
  sequence as values. Index must be the path to a fasta file.'''
  p = subprocess.Popen(args=["bowtie-inspect","-a", "0", index], stdout=subprocess.PIPE, stderr=subprocess.STDOUT) # bowtie-inspect outputs sequences on single lines
  outputlines = p.stdout.readlines()
  p.wait()
  item_dic = {}
  for line in outputlines:
    if (line[0] == ">"):
      try:
        item_dic[current_item] = "".join(stringlist) # to dump the sequence of the previous item - try because of the keyerror of the first item
      except: pass
      current_item = line[1:].rstrip().split()[0] #take the first word before space because bowtie splits headers !
      item_dic[current_item] = ""
      stringlist=[]
    else:
      stringlist.append(line.rstrip() )
  item_dic[current_item] = "".join(stringlist) # for the last item
  return item_dic

def get_fasta_headers (index):
  p = subprocess.Popen(args=["bowtie-inspect","-n", index], stdout=subprocess.PIPE, stderr=subprocess.STDOUT) # bowtie-inspect outputs sequences on single lines
  outputlines = p.stdout.readlines()
  p.wait()
  item_dic = {}
  for line in outputlines:
    header = line.rstrip().split()[0] #take the first word before space because bowtie splits headers !
    item_dic[header] = 1
  return item_dic


def get_file_sample (file, numberoflines):
  '''import random to use this function'''
  F=open(file)
  fullfile = F.read().splitlines()
  F.close()
  if len(fullfile) < numberoflines:
    return "sample size exceeds file size"
  return random.sample(fullfile, numberoflines)

def get_fasta_from_history (file):
  F = open (file, "r")
  item_dic = {}
  for line in F:
    if (line[0] == ">"):
      try:
        item_dic[current_item] = "".join(stringlist) # to dump the sequence of the previous item - try because of the keyerror of the first item
      except: pass
      current_item = line[1:-1].split()[0] #take the first word before space because bowtie splits headers !
      item_dic[current_item] = ""
      stringlist=[]
    else:
      stringlist.append(line[:-1])
  item_dic[current_item] = "".join(stringlist) # for the last item
  return item_dic

def antipara (sequence):
    antidict = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    revseq = sequence[::-1]
    return "".join([antidict[i] for i in revseq])

def RNAtranslate (sequence):
    return "".join([i if i in "AGCN" else "U" for i in sequence])
def DNAtranslate (sequence):
    return "".join([i if i in "AGCN" else "T" for i in sequence])

def RNAfold (sequence_list):
  thestring= "\n".join(sequence_list)
  p = subprocess.Popen(args=["RNAfold","--noPS"], stdin= subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  output=p.communicate(thestring)[0]
  p.wait()
  output=output.split("\n")
  if not output[-1]: output = output[:-1] # nasty patch to remove last empty line
  buffer=[]
  for line in output:
    if line[0] in ["N","A","T","U","G","C"]:
      buffer.append(DNAtranslate(line))
    if line[0] in ["(",".",")"]:
      fields=line.split("(")
      energy= fields[-1]
      energy = energy[:-1] # remove the ) parenthesis
      energy=float(energy)
      buffer.append(str(energy))
  return dict(zip(buffer[::2], buffer[1::2]))

def extractsubinstance (start, end, instance):
  ''' Testing whether this can be an function external to the class to save memory'''
  subinstance = SmRNAwindow (instance.gene, instance.sequence[start-1:end], start)
  subinstance.gene = "%s %s %s" % (subinstance.gene, subinstance.windowoffset, subinstance.windowoffset + subinstance.size - 1)
  upcoordinate = [i for i in range(start,end+1) if instance.readDict[i] ]
  downcoordinate = [-i for i in range(start,end+1) if instance.readDict[-i] ]
  for i in upcoordinate:
    subinstance.readDict[i]=instance.readDict[i]
  for i in downcoordinate:
    subinstance.readDict[i]=instance.readDict[i]
  return subinstance

class HandleSmRNAwindows:
  def __init__(self, alignmentFile="~", alignmentFileFormat="tabular", genomeRefFile="~", genomeRefFormat="bowtieIndex", biosample="undetermined"):
    self.biosample = biosample
    self.alignmentFile = alignmentFile
    self.alignmentFileFormat = alignmentFileFormat # can be "tabular" or "sam"
    self.genomeRefFile = genomeRefFile
    self.genomeRefFormat = genomeRefFormat # can be "bowtieIndex" or "fastaSource"
    self.alignedReads = 0
    self.instanceDict = {}
    if genomeRefFormat == "bowtieIndex":
      self.itemDict = get_fasta (genomeRefFile)
    elif genomeRefFormat == "fastaSource":
      self.itemDict = get_fasta_from_history (genomeRefFile)
    for item in self.itemDict:
      self.instanceDict[item] = SmRNAwindow(item, sequence=self.itemDict[item], windowoffset=1, biosample=self.biosample) # create as many instances as there is items
    self.readfile()

  def readfile (self) :
    if self.alignmentFileFormat == "tabular":
      F = open (self.alignmentFile, "r")
      for line in F:
        fields = line.split()
        polarity = fields[1]
        gene = fields[2]
        offset = int(fields[3])
        size = len (fields[4])
        self.instanceDict[gene].addread (polarity, offset+1, size) # to correct to 1-based coordinates of SmRNAwindow
        self.alignedReads += 1
      F.close()
    elif self.alignmentFileFormat == "sam":
      F = open (self.alignmentFile, "r")
      dict = {"0":"+", "16":"-"}
      for line in F:
        if line[0]=='@':
            continue
        fields = line.split()
        if fields[2] == "*": continue
        polarity = dict[fields[1]]
        gene = fields[2]
        offset = int(fields[3])
        size = len (fields[9])
        self.instanceDict[gene].addread (polarity, offset, size) # sam format is already 1-based coordinates
        self.alignedReads += 1
      F.close()
    elif self.alignmentFileFormat == "bam":
      import pysam
      samfile = pysam.Samfile(self.alignmentFile)
      for read in samfile:
        if read.tid == -1:
          continue # filter out unaligned reads
        if read.is_reverse:
          polarity="-"
        else:
          polarity="+"
        gene = samfile.getrname(read.tid)
        offset = read.pos
        size = read.qlen
        self.instanceDict[gene].addread (polarity, offset+1, size) # pysam converts coordinates to 0-based (https://media.readthedocs.org/pdf/pysam/latest/pysam.pdf)
        self.alignedReads += 1
      return

  def CountFeatures (self, GFF3="path/to/file"):
    featureDict = defaultdict(int)
    F  = open (GFF3, "r")
    for line in F:
      if line[0] ==  "#": continue
      fields = line[:-1].split()
      chrom, feature, leftcoord, rightcoord, polarity = fields[0], fields[2], fields[3], fields[4], fields[6]
      featureDict[feature] += self.instanceDict[chrom].readcount(upstream_coord=int(leftcoord), downstream_coord=int(rightcoord), polarity="both", method="destructive")
    F.close()
    return featureDict

class SmRNAwindow:

  def __init__(self, gene, sequence="ATGC", windowoffset=1, biosample="Undetermined"):
    self.biosample = biosample
    self.sequence = sequence
    self.gene = gene
    self.windowoffset = windowoffset
    self.size = len(sequence)
    self.readDict = defaultdict(list) # with a {+/-offset:[size1, size2, ...], ...}
    self.matchedreadsUp = 0
    self.matchedreadsDown = 0
    
  def addread (self, polarity, offset, size):
    '''ATTENTION ATTENTION ATTENTION'''
    ''' We removed the conversion from 0 to 1 based offset, as we do this now during readparsing.'''
    if polarity == "+":
      self.readDict[offset].append(size)
      self.matchedreadsUp += 1
    else:
      self.readDict[-(offset + size -1)].append(size)
      self.matchedreadsDown += 1
    return

  def barycenter (self, upstream_coord=None, downstream_coord=None):
    '''refactored 24-12-2013 to save memory and introduce offset filtering see readcount method for further discussion on that
    In this version, attempt to replace the dictionary structure by a list of tupple to save memory too'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    window_size = downstream_coord - upstream_coord +1
    def weigthAverage (TuppleList):
      weightSum = 0
      PonderWeightSum = 0
      for tuple in TuppleList:
        PonderWeightSum += tuple[0] * tuple[1]
        weightSum += tuple[1]
      if weightSum > 0:
        return PonderWeightSum / float(weightSum)
      else:
        return 0
    forwardTuppleList = [(k, len(self.readDict[k])) for k in self.readDict.keys() if (k > 0 and abs(k) >= upstream_coord and abs(k) <= downstream_coord)] # both forward and in the proper offset window
    reverseTuppleList = [(-k, len(self.readDict[k])) for k in self.readDict.keys() if (k < 0 and abs(k) >= upstream_coord and abs(k) <= downstream_coord)] # both reverse and in the proper offset window
    Fbarycenter = (weigthAverage (forwardTuppleList) - upstream_coord) / window_size
    Rbarycenter = (weigthAverage (reverseTuppleList) - upstream_coord) / window_size
    return Fbarycenter, Rbarycenter

  def correlation_mapper (self, reference, window_size):
    '''to map correlation with a sliding window 26-2-2013'''
    if window_size > self.size:
      return []
    F=open(reference, "r")
    reference_forward = []
    reference_reverse = []
    for line in F:
      fields=line.split()
      reference_forward.append(int(float(fields[1])))
      reference_reverse.append(int(float(fields[2])))
    F.close()
    local_object_forward=[]
    local_object_reverse=[]
    ## Dict to list for the local object
    for i in range(1, self.size+1):
      local_object_forward.append(len(self.readDict[i]))
      local_object_reverse.append(len(self.readDict[-i]))
    ## start compiling results by slides
    results=[]
    for coordinate in range(self.size - window_size):
      local_forward=local_object_forward[coordinate:coordinate + window_size]
      local_reverse=local_object_reverse[coordinate:coordinate + window_size]
      if sum(local_forward) == 0 or sum(local_reverse) == 0: 
        continue
      try:
        reference_to_local_cor_forward = stats.spearmanr(local_forward, reference_forward)
        reference_to_local_cor_reverse = stats.spearmanr(local_reverse, reference_reverse)
        if (reference_to_local_cor_forward[0] > 0.2 or  reference_to_local_cor_reverse[0]>0.2):
          results.append([coordinate+1, reference_to_local_cor_forward[0], reference_to_local_cor_reverse[0]])
      except:
        pass
    return results

  def readcount (self, size_inf=0, size_sup=1000, upstream_coord=None, downstream_coord=None, polarity="both", method="conservative"):
    '''refactored 24-12-2013 to save memory and introduce offset filtering
    take a look at the defaut parameters that cannot be defined relatively to the instance are they are defined before instanciation
    the trick is to pass None and then test
    polarity parameter can take "both", "forward" or "reverse" as value'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    if upstream_coord == 1 and downstream_coord == self.windowoffset+self.size-1 and polarity == "both":
      return self.matchedreadsUp +  self.matchedreadsDown
    if upstream_coord == 1 and downstream_coord == self.windowoffset+self.size-1 and polarity == "forward":
      return self.matchedreadsUp    
    if upstream_coord == 1 and downstream_coord == self.windowoffset+self.size-1 and polarity == "reverse":
      return self.matchedreadsDown    
    n=0
    if polarity == "both":
      for offset in xrange(upstream_coord, downstream_coord+1):
        if self.readDict.has_key(offset):
          for read in self.readDict[offset]:
            if (read>=size_inf and read<= size_sup):
              n += 1
          if method != "conservative":
            del self.readDict[offset] ## Carefull ! precludes re-use on the self.readDict dictionary !!!!!! TEST
        if self.readDict.has_key(-offset):
          for read in self.readDict[-offset]:
            if (read>=size_inf and read<= size_sup):
              n += 1
          if method != "conservative":
            del self.readDict[-offset]
      return n
    elif polarity == "forward":
      for offset in xrange(upstream_coord, downstream_coord+1):
        if self.readDict.has_key(offset):
          for read in self.readDict[offset]:
            if (read>=size_inf and read<= size_sup):
              n += 1
      return n
    elif polarity == "reverse":
      for offset in xrange(upstream_coord, downstream_coord+1):
        if self.readDict.has_key(-offset):
          for read in self.readDict[-offset]:
            if (read>=size_inf and read<= size_sup):
              n += 1
      return n

  def readsizes (self):
    '''return a dictionary of number of reads by size (the keys)'''
    dicsize = {}
    for offset in self.readDict:
      for size in self.readDict[offset]:
        dicsize[size] = dicsize.get(size, 0) + 1
    return dicsize
    
  def statsizes (self, upstream_coord=None, downstream_coord=None):
    ''' migration to memory saving by specifying possible subcoordinates
    see the readcount method for further discussion'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    L = []
    for offset in self.readDict:
      if (abs(offset) < upstream_coord or abs(offset) > downstream_coord): continue
      for size in self.readDict[offset]:
        L.append(size)
    meansize = mean(L)
    stdv = std(L)
    mediansize = median(L)         
    return meansize, mediansize, stdv

  def foldEnergy (self, upstream_coord=None, downstream_coord=None):
    ''' migration to memory saving by specifying possible subcoordinates
    see the readcount method for further discussion'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    Energy = RNAfold ([self.sequence[upstream_coord-1:downstream_coord] ])
    return float(Energy[self.sequence[upstream_coord-1:downstream_coord]])

  def Ufreq (self, size_scope, upstream_coord=None, downstream_coord=None):
    ''' migration to memory saving by specifying possible subcoordinates
    see the readcount method for further discussion. size_scope must be an interable'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    freqDic = {"A":0,"T":0,"G":0,"C":0, "N":0}
    convertDic = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    for offset in self.readDict:
      if (abs(offset) < upstream_coord or abs(offset) > downstream_coord): continue
      for size in self.readDict[offset]:
        if size in size_scope:
          startbase = self.sequence[abs(offset)-self.windowoffset]
          if offset < 0:
            startbase = convertDic[startbase]
          freqDic[startbase] += 1
    base_sum = float ( sum( freqDic.values()) )
    if base_sum == 0:
      return "."
    else:
      return freqDic["T"] / base_sum * 100

  def Ufreq_stranded (self, size_scope, upstream_coord=None, downstream_coord=None):
    ''' migration to memory saving by specifying possible subcoordinates
    see the readcount method for further discussion. size_scope must be an interable
    This method is similar to the Ufreq method but take strandness into account'''
    upstream_coord = upstream_coord or self.windowoffset
    downstream_coord = downstream_coord or self.windowoffset+self.size-1
    freqDic = {"Afor":0,"Tfor":0,"Gfor":0,"Cfor":0, "Nfor":0,"Arev":0,"Trev":0,"Grev":0,"Crev":0, "Nrev":0}
    convertDic = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    for offset in self.readDict:
      if (abs(offset) < upstream_coord or abs(offset) > downstream_coord): continue
      for size in self.readDict[offset]:
        if size in size_scope:
          startbase = self.sequence[abs(offset)-self.windowoffset]
          if offset < 0:
            startbase = convertDic[startbase]
            freqDic[startbase+"rev"] += 1
          else:
            freqDic[startbase+"for"] += 1
    forward_sum = float ( freqDic["Afor"]+freqDic["Tfor"]+freqDic["Gfor"]+freqDic["Cfor"]+freqDic["Nfor"])
    reverse_sum = float ( freqDic["Arev"]+freqDic["Trev"]+freqDic["Grev"]+freqDic["Crev"]+freqDic["Nrev"])
    if forward_sum == 0 and reverse_sum == 0:
      return ". | ."
    elif reverse_sum == 0:
      return "%s | ." % (freqDic["Tfor"] / forward_sum * 100)
    elif forward_sum == 0:
      return ". | %s" % (freqDic["Trev"] / reverse_sum * 100)
    else:
      return "%s | %s" % (freqDic["Tfor"] / forward_sum * 100, freqDic["Trev"] / reverse_sum * 100)

    
  def readplot (self):
    readmap = {}
    for offset in self.readDict.keys():
      readmap[abs(offset)] = ( len(self.readDict[-abs(offset)]) , len(self.readDict[abs(offset)]) )
    mylist = []
    for offset in sorted(readmap):
      if readmap[offset][1] != 0:
        mylist.append("%s\t%s\t%s\t%s" % (self.gene, offset, readmap[offset][1], "F") )
      if readmap[offset][0] != 0:
        mylist.append("%s\t%s\t%s\t%s" % (self.gene, offset, -readmap[offset][0], "R") )
    return mylist

  def readcoverage (self, upstream_coord=None, downstream_coord=None, windowName=None):
    '''This method has not been tested yet 15-11-2013'''
    upstream_coord = upstream_coord or 1
    downstream_coord = downstream_coord or self.size
    windowName = windowName or "%s_%s_%s" % (self.gene, upstream_coord, downstream_coord)
    forORrev_coverage = dict ([(i,0) for i in xrange(1, downstream_coord-upstream_coord+1)])
    totalforward = self.readcount(upstream_coord=upstream_coord, downstream_coord=downstream_coord, polarity="forward")
    totalreverse = self.readcount(upstream_coord=upstream_coord, downstream_coord=downstream_coord, polarity="reverse")
    if totalforward > totalreverse:
      majorcoverage = "forward"
      for offset in self.readDict.keys():
        if (offset > 0) and ((offset-upstream_coord+1) in forORrev_coverage.keys() ):
          for read in self.readDict[offset]:
            for i in xrange(read):
              try:
                forORrev_coverage[offset-upstream_coord+1+i] += 1
              except KeyError:
                continue # a sense read may span over the downstream limit
    else:
      majorcoverage = "reverse"
      for offset in self.readDict.keys():
        if (offset < 0) and (-offset-upstream_coord+1 in forORrev_coverage.keys() ):
          for read in self.readDict[offset]:
            for i in xrange(read):
              try:
                forORrev_coverage[-offset-upstream_coord-i] += 1 ## positive coordinates in the instance, with + for forward coverage and - for reverse coverage
              except KeyError:
                continue # an antisense read may span over the upstream limit
    output_list = []
    maximum = max (forORrev_coverage.values()) or 1
    for n in sorted (forORrev_coverage):
      output_list.append("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.biosample, windowName, n, float(n)/(downstream_coord-upstream_coord+1), forORrev_coverage[n], float(forORrev_coverage[n])/maximum, majorcoverage))
    return "\n".join(output_list)

          
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

  def phasing (self, size_range, scope):
    ''' to calculate autocorelation like signal - scope must be an python iterable'''
    read_table = {}
    total_read_number = 0
    general_frequency_table = dict ([(i, 0) for i in scope])
    ## read input filtering
    for offset in self.readDict:
      for size in self.readDict[offset]:
        if size in size_range:
          read_table[offset] = read_table.get(offset, 0) + 1
          total_read_number += 1
    ## per offset read phasing computing
    for offset in read_table:
      frequency_table = dict ([(i, 0) for i in scope]) # local frequency table
      number_of_targets = 0
      for i in scope:
        if offset > 0:
          frequency_table[i] += read_table[offset] *  read_table.get(offset + i, 0)
          number_of_targets += read_table.get(offset + i, 0)
        else:
          frequency_table[i] += read_table[offset] *  read_table.get(offset - i, 0)
          number_of_targets += read_table.get(offset - i, 0)
    ## inclusion of local frequency table in the general frequency table (all offsets average)
      for i in scope:
        try:
          general_frequency_table[i] += (1. / number_of_targets / total_read_number) * frequency_table[i]
        except ZeroDivisionError :
          continue
    return general_frequency_table



  def z_signature (self, minquery, maxquery, mintarget, maxtarget, scope):
    '''Must do: from numpy import mean, std, to use this method; scope must be a python iterable and defines the relative offsets to compute'''
    frequency_table = self.signature (minquery, maxquery, mintarget, maxtarget, scope)
    z_table = {}
    frequency_list = [frequency_table[i] for i in sorted (frequency_table)]
    if std(frequency_list):
      meanlist = mean(frequency_list)
      stdlist = std(frequency_list)
      z_list = [(i-meanlist)/stdlist for i in frequency_list]
      return dict (zip (sorted(frequency_table), z_list) ) 
    else:
      return dict (zip (sorted(frequency_table), [0 for i in frequency_table]) )

  def percent_signature (self, minquery, maxquery, mintarget, maxtarget, scope):
    frequency_table = self.signature (minquery, maxquery, mintarget, maxtarget, scope)
    total = float(sum ([self.readsizes().get(i,0) for i in set(range(minquery,maxquery)+range(mintarget,maxtarget))]) )
    if total == 0:
      return dict( [(i,0) for i in scope])
    return dict( [(i, frequency_table[i]/total*100) for i in scope])

  def pairer (self, overlap, minquery, maxquery, mintarget, maxtarget):
    queryhash = defaultdict(list)
    targethash = defaultdict(list)
    query_range = range (int(minquery), int(maxquery)+1)
    target_range = range (int(mintarget), int(maxtarget)+1)
    paired_sequences = []
    for offset in self.readDict: # selection of data
      for size in self.readDict[offset]:
        if size in query_range:
          queryhash[offset].append(size)
        if size in target_range:
          targethash[offset].append(size)
    for offset in queryhash:
      if offset >= 0: matched_offset = -offset - overlap + 1
      else: matched_offset = -offset - overlap + 1
      if targethash[matched_offset]:
        paired = min ( len(queryhash[offset]), len(targethash[matched_offset]) )
        if offset >= 0:
          for i in range (paired):
            paired_sequences.append("+%s" % RNAtranslate ( self.sequence[offset:offset+queryhash[offset][i]]) )
            paired_sequences.append("-%s" % RNAtranslate (antipara (self.sequence[-matched_offset-targethash[matched_offset][i]+1:-matched_offset+1]) ) )
        if offset < 0:
          for i in range (paired):
            paired_sequences.append("-%s" % RNAtranslate (antipara (self.sequence[-offset-queryhash[offset][i]+1:-offset+1]) ) )
            paired_sequences.append("+%s" % RNAtranslate (self.sequence[matched_offset:matched_offset+targethash[matched_offset][i]] ) )
    return paired_sequences

  def pairable (self, overlap, minquery, maxquery, mintarget, maxtarget):
    queryhash = defaultdict(list)
    targethash = defaultdict(list)
    query_range = range (int(minquery), int(maxquery)+1)
    target_range = range (int(mintarget), int(maxtarget)+1)
    paired_sequences = []

    for offset in self.readDict: # selection of data
      for size in self.readDict[offset]:
        if size in query_range:
          queryhash[offset].append(size)
        if size in target_range:
          targethash[offset].append(size)

    for offset in queryhash:
      matched_offset = -offset - overlap + 1
      if targethash[matched_offset]:
        if offset >= 0:
          for i in queryhash[offset]:
            paired_sequences.append("+%s" % RNAtranslate (self.sequence[offset:offset+i]) )
          for i in targethash[matched_offset]:
            paired_sequences.append( "-%s" % RNAtranslate (antipara (self.sequence[-matched_offset-i+1:-matched_offset+1]) ) )
        if offset < 0:
          for i in queryhash[offset]:
            paired_sequences.append("-%s" %  RNAtranslate (antipara (self.sequence[-offset-i+1:-offset+1]) ) )
          for i in targethash[matched_offset]:
            paired_sequences.append("+%s" %  RNAtranslate (self.sequence[matched_offset:matched_offset+i] ) )
    return paired_sequences

  def newpairable_bowtie (self, overlap, minquery, maxquery, mintarget, maxtarget):
    ''' revision of pairable on 3-12-2012, with focus on the offset shift problem (bowtie is 1-based cooordinates whereas python strings are 0-based coordinates'''
    queryhash = defaultdict(list)
    targethash = defaultdict(list)
    query_range = range (int(minquery), int(maxquery)+1)
    target_range = range (int(mintarget), int(maxtarget)+1)
    bowtie_output = []

    for offset in self.readDict: # selection of data
      for size in self.readDict[offset]:
        if size in query_range:
          queryhash[offset].append(size)
        if size in target_range:
          targethash[offset].append(size)
    counter = 0
    for offset in queryhash:
      matched_offset = -offset - overlap + 1
      if targethash[matched_offset]:
        if offset >= 0:
          for i in queryhash[offset]:
            counter += 1
            bowtie_output.append("%s\t%s\t%s\t%s\t%s" % (counter, "+", self.gene, offset-1, self.sequence[offset-1:offset-1+i]) ) # attention a la base 1-0 de l'offset 
        if offset < 0:
          for i in queryhash[offset]:
            counter += 1
            bowtie_output.append("%s\t%s\t%s\t%s\t%s" % (counter, "-", self.gene, -offset-i, self.sequence[-offset-i:-offset])) # attention a la base 1-0 de l'offset
    return bowtie_output


def __main__(bowtie_index_path, bowtie_output_path):
  sequenceDic = get_fasta (bowtie_index_path)
  objDic = {}
  F = open (bowtie_output_path, "r") # F is the bowtie output taken as input
  for line in F:
    fields = line.split()
    polarity = fields[1]
    gene = fields[2]
    offset = int(fields[3])
    size = len (fields[4])
    try:
      objDic[gene].addread (polarity, offset, size)
    except KeyError:
      objDic[gene] = SmRNAwindow(gene, sequenceDic[gene])
      objDic[gene].addread (polarity, offset, size)
  F.close()
  for gene in objDic:
    print gene, objDic[gene].pairer(19,19,23,19,23)

if __name__ == "__main__" : __main__(sys.argv[1], sys.argv[2]) 
