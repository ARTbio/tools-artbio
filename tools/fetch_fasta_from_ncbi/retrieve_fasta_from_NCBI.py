#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
From a taxonomy ID retrieves all the nucleotide sequences
It returns a multiFASTA nuc/prot file

Entrez Database  UID common name  E-utility Database Name
Nucleotide       GI number        nuccore
Protein          GI number        protein

Retrieve strategy:

esearch to get total number of UIDs (count)
esearch to get UIDs in batches
loop untile end of UIDs list:
  epost to put a batch of UIDs in the history server
  efetch to retrieve info from previous post

retmax of efetch is 1/10 of declared value from NCBI

queries are 1 sec delayed, to satisfy NCBI guidelines (more than what they request)


"""
import sys
import logging
import optparse
import time
import urllib
import urllib2
import httplib
import re


class QueryException(Exception):
    pass


class Eutils:

    def __init__(self, options, logger):
        """
        Initialize retrieval parameters 
        """
        self.logger = logger
        self.base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.query_string = options.query_string
        self.dbname = options.dbname
        if options.outname:
            self.outname = options.outname
        else:
            self.outname = 'NCBI_download' + '.' + self.dbname + '.fasta'
        self.ids = []
        self.retmax_esearch = 100000
        self.retmax_efetch = 500
        self.count = 0
        self.webenv = ""
        self.query_key = ""

    def retrieve(self):
        """
        Retrieve the fasta sequences corresponding to the query
        """
        self.get_count_value()

        # If no UIDs are found exit script
        if self.count > 0:
            self.get_uids_list()
            try:
                self.get_sequences()
            except QueryException as e:
                self.logger.error("Exiting script.")
                raise e
        else:
            self.logger.error("No UIDs were found. Exiting script.")
            raise Exception("")

    def get_count_value(self):
        """
        just to retrieve Count (number of UIDs)
        Total number of UIDs from the retrieved set to be shown in the XML
        output (default=20). By default, ESearch only includes the first 20
        UIDs retrieved in the XML output. If usehistory is set to 'y',
        the remainder of the retrieved set will be stored on the History server;

        http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        """
        self.logger.info("retrieving data from %s" % self.base)
        self.logger.info("for Query: %s and database: %s" %
                         (self.query_string, self.dbname))
        querylog = self.esearch(self.dbname, self.query_string, '', '', "count")
        self.logger.debug("Query response:")
        for line in querylog:
            self.logger.debug(line.rstrip())
            if '</Count>' in line:
                self.count = int(line[line.find('<Count>')+len('<Count>') : line.find('</Count>')])
        self.logger.info("Found %d UIDs" % self.count)

    def get_uids_list(self):
        """
        Increasing retmax allows more of the retrieved UIDs to be included in the XML output,
        up to a maximum of 100,000 records.
        from http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
        """
        retmax = self.retmax_esearch
        if (self.count > retmax):
            num_batches = (self.count / retmax) + 1
        else:
            num_batches = 1
        self.logger.info("Batch size for esearch action: %d UIDs" % retmax)
        self.logger.info("Number of batches for esearch action: %d " % num_batches)
        for n in range(num_batches):
            querylog = self.esearch(self.dbname, self.query_string, n*retmax, retmax, '')
            for line in querylog:
                if '<Id>' in line and '</Id>' in line:
                    uid = (line[line.find('<Id>')+len('<Id>') : line.find('</Id>')])
                    self.ids.append(uid)
            self.logger.info("Retrieved %d UIDs" % len(self.ids))

    def esearch(self, db, term, retstart, retmax, rettype):
        url = self.base + "esearch.fcgi"
        self.logger.debug("url: %s" % url)
        values = {'db': db,
                  'term': term,
                  'rettype': rettype,
                  'retstart': retstart,
                  'retmax': retmax}
        data = urllib.urlencode(values)
        self.logger.debug("data: %s" % str(data))
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        querylog = response.readlines()
        response.close()
        time.sleep(1)
        return querylog

    def epost(self, db, ids):
        url = self.base + "epost.fcgi"
        self.logger.debug("url_epost: %s" % url)
        values = {'db': db,
                  'id': ids}
        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        serverResponse = False
        nb_trials = 0
        while not serverResponse:
            nb_trials += 1
            try:
                self.logger.debug("Try number %s for opening and readin URL %s" % ( nb_trials, url+data ))
                response = urllib2.urlopen(req)
                querylog = response.readlines()
                response.close()
                serverResponse = True
            except urllib2.HTTPError as e:
                self.logger.info("urlopen error:%s, %s" % (e.code, e.read() ) )
                self.logger.info("Retrying in 1 sec")
                serverResponse = False
                time.sleep(1)
            except urllib2.URLError as e:
                self.logger.info("urlopen error: Failed to reach a server")
                self.logger.info("Reason :%s" % ( e.reason ) )
                self.logger.info("Retrying in 1 sec")
                serverResponse = False
                time.sleep(1)
            except httplib.IncompleteRead as e:
                self.logger.info("IncompleteRead error:  %s" % ( e.partial ) )
                self.logger.info("Retrying in 1 sec")
                serverResponse = False
                time.sleep(1)
        self.logger.debug("query response:")
        for line in querylog:
            self.logger.debug(line.rstrip())
            if '</QueryKey>' in line:
                self.query_key = str(line[line.find('<QueryKey>')+len('<QueryKey>'):line.find('</QueryKey>')])
            if '</WebEnv>' in line:
                self.webenv = str(line[line.find('<WebEnv>')+len('<WebEnv>'):line.find('</WebEnv>')])
            self.logger.debug("*** epost action ***")
            self.logger.debug("query_key: %s" % self.query_key)
            self.logger.debug("webenv: %s" % self.webenv)
        time.sleep(1)

    def efetch(self, db, query_key, webenv):
        url = self.base + "efetch.fcgi"
        self.logger.debug("url_efetch: %s" % url)
        values = {'db': db,
                  'query_key': query_key,
                  'webenv': webenv,
                  'rettype': "fasta",
                  'retmode': "text"}
        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        self.logger.debug("data: %s" % str(data))
        serverTransaction = False
        counter = 0
        response_code = 0
        while not serverTransaction:
            counter += 1
            self.logger.info("Server Transaction Trial:  %s" % ( counter ) )
            try:
                response = urllib2.urlopen(req)
                response_code = response.getcode()
                fasta = response.read()
                response.close()
                if ( (response_code != 200) or ("Resource temporarily unavailable" in fasta)
                    or ("Error" in fasta) or (not fasta.startswith(">") ) ):
                    serverTransaction = False
                    if ( response_code != 200 ):
                        self.logger.info("urlopen error: Response code is not 200")
                    elif ( "Resource temporarily unavailable" in fasta ):
                        self.logger.info("Ressource temporarily unavailable")
                    elif ( "Error" in fasta ):
                        self.logger.info("Error in fasta")
                    else:
                        self.logger.info("Fasta doesn't start with '>'")
                else:
                    serverTransaction = True
            except urllib2.HTTPError as e:
                serverTransaction = False
                self.logger.info("urlopen error:%s, %s" % (e.code, e.read() ) )
            except urllib2.URLError as e:
                serverTransaction = False
                self.logger.info("urlopen error: Failed to reach a server")
                self.logger.info("Reason :%s" % ( e.reason ) )
            except httplib.IncompleteRead as e:
                serverTransaction = False
                self.logger.info("IncompleteRead error:  %s" % ( e.partial ) )
            if (counter > 500):
                serverTransaction = True
        if (counter > 500):
            raise QueryException({"message":"500 Server Transaction Trials attempted for this batch. Aborting."})
        fasta = self.sanitiser(self.dbname, fasta) 
        time.sleep(0.1)
        return fasta
        
    def sanitiser(self, db, fastaseq):
        if db not in "nuccore protein" : return fastaseq
        regex = re.compile(r"[ACDEFGHIKLMNPQRSTVWYBZ]{49,}")
        sane_seqlist = []
        seqlist = fastaseq.split("\n\n")
        for seq in seqlist[:-1]:
            fastalines = seq.split("\n")
            if len(fastalines) < 2:
                self.logger.info("Empty sequence for %s" % ("|".join(fastalines[0].split("|")[:4]) ) )
                self.logger.info("%s download is skipped" % ("|".join(fastalines[0].split("|")[:4]) ) )
                continue
            if db == "nuccore":
                badnuc = 0
                for nucleotide in fastalines[1]:
                    if nucleotide not in "ATGC":
                        badnuc += 1
                if float(badnuc)/len(fastalines[1]) > 0.4:
                    self.logger.info("%s ambiguous nucleotides in %s or download interrupted at this offset | %s" % ( float(badnuc)/len(fastalines[1]), "|".join(fastalines[0].split("|")[:4]), fastalines[1]) )
                    self.logger.info("%s download is skipped" % (fastalines[0].split("|")[:4]) )
                    continue
                fastalines[0] = fastalines[0].replace(" ","_")[:100] # remove spaces and trim the header to 100 chars
                cleanseq = "\n".join(fastalines)
                sane_seqlist.append(cleanseq)
            elif db == "protein":
                fastalines[0] = fastalines[0][0:100]
                fastalines[0] = fastalines[0].replace(" ", "_")
                fastalines[0] = fastalines[0].replace("[", "_")
                fastalines[0] = fastalines[0].replace("]", "_")
                fastalines[0] = fastalines[0].replace("=", "_")
                fastalines[0] = fastalines[0].rstrip("_") # because blast makedb doesn't like it 
                fastalines[0] = re.sub(regex, "_", fastalines[0])
                cleanseq = "\n".join(fastalines)
                sane_seqlist.append(cleanseq)
        self.logger.info("clean sequences appended: %d" % (len(sane_seqlist) ) )
        return "\n".join(sane_seqlist)

    def get_sequences(self):
        """
        Total number of records from the input set to be retrieved, up to a maximum
        of 10,000. Optionally, for a large set the value of retstart can be iterated
        while holding retmax constant, thereby downloading the entire set in batches
        of size retmax.
        
        http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        
        """
        batch_size = self.retmax_efetch
        count = self.count
        uids_list = self.ids
        self.logger.info("Batch size for efetch action: %d" % batch_size)
        self.logger.info("Number of batches for efetch action: %d" % ((count / batch_size) + 1))
        with open(self.outname, 'w') as out:
            for start in range(0, count, batch_size):
                end = min(count, start+batch_size)
                batch = uids_list[start:end]
                if self.epost(self.dbname, ",".join(batch)) != -1:
                    mfasta = ''
                    while not mfasta:
                        self.logger.info("retrieving batch %d" % ((start / batch_size) + 1))
                        try:
                            mfasta = self.efetch(self.dbname, self.query_key, self.webenv)
                            out.write(mfasta + '\n')
                        except QueryException as e:
                            self.logger.error("%s" % e.message)
                            raise e


LOG_FORMAT = '%(asctime)s|%(levelname)-8s|%(message)s'
LOG_DATEFMT = '%Y-%m-%d %H:%M:%S'
LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


def __main__():
    """ main function """
    parser = optparse.OptionParser(description='Retrieve data from NCBI')
    parser.add_option('-i', dest='query_string', help='NCBI Query String')
    parser.add_option('-o', dest='outname', help='output file name')
    parser.add_option('-l', '--logfile', help='log file (default=stderr)')
    parser.add_option('--loglevel', choices=LOG_LEVELS, default='INFO', help='logging level (default: INFO)')
    parser.add_option('-d', dest='dbname', help='database type')
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')
    
    log_level = getattr(logging, options.loglevel)
    kwargs = {'format': LOG_FORMAT,
              'datefmt': LOG_DATEFMT,
              'level': log_level}
    if options.logfile:
        kwargs['filename'] = options.logfile
    logging.basicConfig(**kwargs)
    logger = logging.getLogger('data_from_NCBI')
    
    E = Eutils(options, logger)
    try:
        E.retrieve()
    except Exception as e:
        sys.exit(1)


if __name__ == "__main__":
    __main__()
