#!/usr/bin/python3

import argparse


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--input', action="store", type=str,
        help="varscan vcf file with normal and tumor genotypes (columns 10 and 11)")
    the_parser.add_argument(
        '--output', action="store", type=str,
        help="vcf with computed VAFs")
    args = the_parser.parse_args()
    return args


def main(input, output):
    with open(sys.argv[1], 'r') as f:
        myinput = f.read()
    mylines = myinput.split('\n')
    entete = [i for i in mylines[:-1] if i[0] == '#']
    variant = [i for i in mylines[:-1] if i[0] != '#']
    out = open(output, 'w')
    out.write('\n'.join(entete[:-1]))
    out.write('\n')
    out.write('##FORMAT=<ID=VAF,Number=R,Type=float,Description="Variant \
               Allele Frequency">\n')
    out.write(entete[-1], '\n')
    for i in variant:
        fields = i.split('\t')[9:11]
        af_normal = fields[0].split(':')[3]
        vac_normal = af_normal.split(',')
        af_tumor = fields[1].split(':')[3]
        vac_tumor = af_tumor.split(',')
        vaf_normal = int(vac_normal[1]) / (int(vac_normal[0]) + int(vac_normal[1]))
        vaf_tumor = int(vac_tumor[1]) / (int(vac_tumor[0]) + int(vac_tumor[1]))
        normal_list = fields[0].split(':')
        normal_list.append(f'{vaf_normal:.3f}')
        tumor_list = fields[1].split(':')
        tumor_list.append(f'{vaf_tumor:.3f}')
        normal_string = ':'.join(normal_list)
        tumor_string = ':'.join(tumor_list)
        fields = i.split('\t')[0:9]
        fields[8] += ':VAF'
        fields.append(normal_string)
        fields.append(tumor_string)
        out.write('\t'.join(fields), '\n')
        out.close()


if __name__ == "__main__":
    args = Parser()
    main(args.input, args.output)
