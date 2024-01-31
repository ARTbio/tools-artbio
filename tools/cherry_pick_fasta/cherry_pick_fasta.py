import argparse
from collections import defaultdict


def Parser():
    the_parser = argparse.ArgumentParser(
        description='Cherry pick fasta sequences')
    the_parser.add_argument('--input', action='store', type=str,
                            help='input fasta file')
    the_parser.add_argument('--searchfor', action='store', type=str,
                            help='with, without, or withlist, withoutlist')
    the_parser.add_argument('--mode', action='store', type=str,
                            default='includes', help='exact or includes')
    the_parser.add_argument('--query-string', dest='query_string',
                            action='store', type=str,
                            help='headers containing the string will be \
                                  extracted or excluded as well as the \
                                  corresponding sequence')
    the_parser.add_argument('--query-file', dest='query_file',
                            action='store', type=str,
                            help='headers containing any of the strings \
                                  provided in the text file (1 string per \
                                  line) will be extracted or excluded as well \
                                   as the corresponding sequence')
    the_parser.add_argument('--output', action='store', type=str,
                            help='output fasta file')
    args = the_parser.parse_args()
    return args


def parse_fasta_dict(query, fasta_dict, mode):

    if not isinstance(query, list):
        query = [query]

    def kmers(string, ksize, index):
        if ksize > len(string):
            return
        for i in range(len(string) - ksize + 1):
            kmer = string[i:i+ksize]
            index[kmer].append(string)

    def consult_index(word, index):
        accumulator = []
        print(len(index[word]))
        for title in index[word]:
            accumulator.append(title)
        print(len(accumulator))
        for title in set(accumulator):
            print(title)

    accumulator = []
    if mode == 'includes':
        kmersizes = set([len(word) for word in query])
        index = defaultdict(list)
        for size in kmersizes:
            for header in fasta_dict:
                kmers(header, size, index)
        for keyword in query:
            for header in index[keyword]:
                accumulator.append(header)
        accumulator = set(accumulator)
        res_dict = {k: fasta_dict[k] for k in fasta_dict if k in accumulator}
        return res_dict
    elif mode == 'exact':
        for keyword in query:
            try:
                len(fasta_dict[keyword])
                accumulator.append(keyword)
            except KeyError:
                pass
        accumulator = set(accumulator)
        res_dict = {k: fasta_dict[k] for k in fasta_dict if k in accumulator}
        return res_dict


def complement_fasta_dict(fasta_dict, subfasta_dict):
    fasta_ids = list(fasta_dict.keys())
    subfasta_ids = list(subfasta_dict.keys())
    complement_ids = list(set(fasta_ids) - set(subfasta_ids))
    sub_dict = {k: fasta_dict[k] for k in fasta_dict if k in complement_ids}
    return sub_dict


def getquerylist(file):
    querylist = []
    for line in open(file, 'r'):
        querylist.append(line.rstrip())
    return querylist


def buid_fasta_dict(fasta):
    seq_dict = dict()
    f = open(fasta, 'r')
    content = f.read()
    segmented_content = content.split('>')
    segmented_content = segmented_content[1:]
    for seq in segmented_content:
        sliced_seq = seq.split('\n')
        header = sliced_seq[0]
        sliced_seq = sliced_seq[1:]
        sequence = ''.join(sliced_seq)
        seq_dict[header] = sequence
    return seq_dict


def write_fasta_result(fasta_dict, file):
    line_length = 60
    with open(file, 'w') as f:
        for header in sorted(fasta_dict):
            f.write('>%s\n' % header)
            if len(fasta_dict[header]) <= line_length:
                f.write('%s\n' % fasta_dict[header])
            else:
                for i in range(line_length, len(fasta_dict[header]),
                               line_length):
                    f.write('%s\n' % fasta_dict[header][i-line_length:i])
                f.write('%s\n' % fasta_dict[header][i:])


def __main__():
    ''' main function '''
    args = Parser()
    fasta_dict = buid_fasta_dict(args.input)
    if args.query_string:
        query = args.query_string
    elif args.query_file:
        query = getquerylist(args.query_file)
    if args.searchfor == 'with':
        fasta_result_dict = parse_fasta_dict(query, fasta_dict, args.mode)
    elif args.searchfor == 'without':
        fasta_result_dict = complement_fasta_dict(fasta_dict, parse_fasta_dict(
                                                  query, fasta_dict,
                                                  args.mode))
    write_fasta_result(fasta_result_dict, args.output)


if __name__ == '__main__':
    __main__()
