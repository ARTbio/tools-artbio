#!/usr/bin/env python3
#Author: Leonardo G. Panunzi, 2019

from __future__ import print_function
from lib_repo.scripts import SigProfilerMatrixGeneratorFunc as create_mat
from lib_repo import install as gnm_install
import argparse
import sys

def run_pipe (install_genome,name,genome,files,exome,bed,chrom,plot,tsb,gs):
    if install_genome != None:
        print("Downloading the following reference genome:\t" + str(install_genome), "\n")
        gnm_install.install(install_genome, rsync=False, bash=True)
    else:
        print("Input data:\nPrefix: " + name + "\nReference Genome :" + genome + "\nVCF files:" + files, "\n")
        print("Output results:\t\t" + files + "output/","\n")
        matrices = create_mat.SigProfilerMatrixGeneratorFunc(name,genome,files,exome,bed,chrom,plot,tsb,gs)

parser = argparse.ArgumentParser(
        prog='sigprofiler_wrapper',
        #usage='sigprofiler_wrapper <command> <options>',
        description='Provide the necessary arguments to run Sigprofiler => Mutational Signature from Cancer DNA',
)

parser.add_argument("--install_genome", "-ig", help="Install de novo any of the following reference genomes: 'GRCh37', 'GRCh38', 'mm9' or 'mm10'.")
parser.add_argument("--name", "-n", help="Provide a project name")
parser.add_argument("--genome", "-g", help="Provide a reference genome (ex: GRCh37, GRCh38, mm9 or mm10).")
parser.add_argument("--files", "-f", help="Path where the input vcf files are located.")
parser.add_argument("--exome", "-e", action='store_true', help="Flag to use only the exome or not.")
parser.add_argument("--bed", "-b", help="BED file that contains a list of ranges to be used in generating the matrices.")
parser.add_argument("--chrom", "-c", action='store_true', help="Create the matrices on a per chromosome basis.")
parser.add_argument("--plot", "-p", action='store_true', help="Generate the plots for each context.")
parser.add_argument("--tsb", "-t", action='store_true', help="Performs a transcriptional strand bias test for the 24, 384, and 6144 contexts.")
parser.add_argument("--gs", "-s", action='store_true', help="Performs a gene strand bias test.")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args=parser.parse_args()

install_genome = args.install_genome
name = args.name
genome = args.genome
files = args.files
exome = args.exome
bed = args.bed
chrom = args.chrom
plot = args.plot
tsb = args.tsb
gs = args.gs

run_pipe(install_genome,name,genome,files,exome,bed,chrom,plot,tsb,gs)
