import argparse
import gzip


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--hairpins_path', action="store", type=str,
        help="BASE url. ex: /pub/mirbase/22/")
    the_parser.add_argument(
        '--output', action="store", type=str,
        help="parsed hairpin output in fasta format")
    the_parser.add_argument(
        '--basename', action="store", type=str,
        help="genome basename of the parsed fasta")
    args = the_parser.parse_args()
    return args


def get_fasta_dic(gzipfile):
    '''
    gzipfile value example : 'mirbase/22/hairpin.fa.gz'
    '''
    item_dic = {}
    with gzip.open(gzipfile, 'rb') as f:
        current_item = ''
        stringlist = []
        for line in f:
            line = line.decode('utf-8').strip('\n')
            if (line[0] == ">"):
                # dump the sequence of the previous item
                if current_item and stringlist:
                    item_dic[current_item] = "".join(stringlist)
                # take first word of item '''
                current_item = line[1:].split()[0]
                stringlist = []
            else:
                stringlist.append(line)
        item_dic[current_item] = "".join(stringlist)  # for the last item
    return item_dic


def convert_and_print_hairpins(gzipfile, basename, fasta_output):
    raw_fasta_dict = get_fasta_dic(gzipfile)
    parsed_fasta_dict = {}
    trs = str.maketrans("uU", "tT")
    for head in raw_fasta_dict:
        if basename in head:
            parsed_fasta_dict[head] = raw_fasta_dict[head].translate(trs)
    with open(fasta_output, "w") as output:
        for head in sorted(parsed_fasta_dict):
            output.write('>%s\n%s\n' % (head, parsed_fasta_dict[head]))


def main(hairpins_path, basename, outfile):
    convert_and_print_hairpins(hairpins_path, basename, outfile)


if __name__ == "__main__":
    args = Parser()
    main(args.hairpins_path, args.basename, args.output)
