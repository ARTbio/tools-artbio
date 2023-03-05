#!/usr/bin/env python

import argparse
from collections import defaultdict


def command_parse():
    parser = argparse.ArgumentParser(description='Purge GTF file from genes \
        that are on several chromosomes and list them in a log file')
    parser.add_argument(
        '-i', '--input', dest='input', help='input GTF file', required=True)
    parser.add_argument('-o', '--output', dest='output', help='output file \
        name', default='output.gtf')
    parser.add_argument('-l', '--log', dest='log', help='log of purged \
        genes', default='purged_genes.log')
    args = parser.parse_args()
    return args


def get_genes(gtf_file):
    genes = defaultdict(list)
    with open(gtf_file, 'r') as fh:
        for line in fh:
            if line[0] != '#':
                fields = line[:-1].split("\t")
                chrom = fields[0]
                name_gene = fields[-1].split('gene_id "')[-1].split('"; \
                    transcript_id')[0]
                genes[name_gene].append(chrom)
    return genes


def generate_output(genes, log_file):
    '''
    Search for all genes that are present on several chromosomes. This function
    return a list of these genes in target_genes. It also generate a log tab
    delimited file with one gene per line and with its list of chromosomes
    (coma delimited)
    '''
    output = open(log_file, 'w')
    # output.write('#all genes on several chromosomes' + '\n')
    target_genes = list()
    for name_gene in sorted(genes.keys()):
        genes[name_gene] = set(genes[name_gene])
        if len(genes[name_gene]) > 1:
            target_genes.append(name_gene)
            new_line = '\t'.join([name_gene, ','.join(genes[name_gene])])
            output.write("%s\n" % new_line)
    output.close()
    return target_genes


def purge_gtf(target_genes, gtf_file, output_file):
    '''
    Remove all lines of the gtf file where the gene_id is gene of target_genes
    list.
    '''
    output_gtf = open(output_file, 'w')
    with open(gtf_file, 'r') as gtf_handler:
        for line in gtf_handler:
            fields = line[:-1].split("\t")
            gene_name = fields[-1].split('gene_id "')[-1].split('"; \
                transcript_id')[0]
            if gene_name not in target_genes:
                output_gtf.write(line)
    output_gtf.close()


def __main__():
    args = command_parse()
    genes = get_genes(args.input)
    target_genes = generate_output(genes, args.log)
    purge_gtf(target_genes, args.input, args.output)


if __name__ == "__main__":
    __main__()
