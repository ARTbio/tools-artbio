#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import pysam
#from multiprocessing import Process, Queue


def Parser():
    parser = argparse.ArgumentParser(description='miRNAs counts and coverages')
    parser.add_argument('-a', '--alignment', metavar='FILE', type=str,
                        dest='alignment_file', help='Alignment bam file')
    parser.add_argument('--gff', metavar='FILE', type=str, dest='gff_file',
                        help='GFF3 describing both pre-miRNAs\
                              and mature miRNAs')
    parser.add_argument('-q', '--quality_threshold', type=int,
                        dest='quality_threshold',
                        help='Quality threshold for coverage (default=10)',
                        default=10)
    parser.add_argument('-p', '--pre_mirs', type=str, dest='pre_mirs',
                        help='pre-miRNAs count file path', metavar='FILE')
    parser.add_argument('-m', '--mirs', type=str, dest='mirs',
                        help='mature miRNA count file path', metavar='FILE')
    parser.add_argument('--PDF', metavar='FILE', type=str, dest='pdf_output',
                        help='Output file for the lattice dataframe.')
    args = parser.parse_args()
    if args.pdf_output and not args.pre_mirs:
        parser.error("To get a plotted coverage the -p option is necessary.")
    return args


def get_pre_mir_counts(bamfile):
    """
    Takes a AlignmentFile object and returns a dictionary of counts for reads
    aligning with pre_mirs (as keys)
    """
    count = dict()
    for ref_name in bamfile.references:
        count[ref_name] = bamfile.count(reference=ref_name)
    return count


def get_pre_mir_coverage(bamfile, quality=10):
    """
    Takes a AlignmentFile object and returns a dictionary of lists
    of coverage along the coordinates of pre_mirs (as keys)
    """
    coverage = dict()
    for ref_name, ref_len in zip(bamfile.references, bamfile.lengths):
        coverage[ref_name] = bamfile.count_coverage(reference=ref_name,
                                                    start=0, end=ref_len,
                                                    quality_threshold=quality)
        """ Add the 4 coverage values """
        coverage[ref_name] = np.array([sum(x) for x in
                                       zip(*coverage[ref_name])])
    return coverage


def get_mir_counts(bamfile, gff_file):
    """
    Takes a AlignmentFile and a gff file and computes for
    each 'miRNA' region of the gff the number of reads that hit it
    returns a dict[mir_name] = count
    """
    counts = dict()
    for line in open(gff_file, 'r'):
        if line[0] != '#' and line != "":
            gff_fields = line[:-1].split("\t")
            if gff_fields[2] == 'miRNA':
                mir_name = gff_fields[0]
                premir_name = gff_fields[8].split('=')[-1]
                mir_start = int(gff_fields[3])
                mir_end = int(gff_fields[4])
                # GFF is 1-based, pysam is 0-based.
                counts[mir_name] = bamfile.count(reference=premir_name,
                                                 start=mir_start-1,
                                                 end=mir_end-1)
    return counts


def write_dataframe_coverage(countdict, outfile):
    """
    Takes a dict[pre_mir reference name] = [coverage list]
    and writes a dataframe with columns:
    <gene_type name>, offset, normoffset, counts and normcounts
    in the outfile
    """
    F = open(outfile, 'w')
    F.write('Mir_hairpin\tOffset\tNorm_offset\tCount\tNorm_count\n')
    for ref in sorted(countdict):
        """
        For each reference name in mirs,
        write the coverage of each of its positions
        """
        maximum = max(countdict[ref])
        reference_length = len(countdict[ref])
        for pos, c in enumerate(countdict[ref]):
            """ Compute and write value for each reference position"""
            F.write('%s\t%s\t%s\t%s\t%s\n' % (ref, str(pos + 1),
                    str(float(pos+1)/reference_length), str(float(c)),
                    str(float(c)/maximum) if maximum != 0 else '0'))
    F.close()

def plot_coverage_par(countdict, outfile):
    """
    Tkes a dict[pre_mir reference name] = [coverage list]
    and plots the coverage of each pre_mir reference name
    """
    pdf = PdfPages('multipage.pdf')
    xax = 0
    yax = 0
    process_list = list()
    nb_images = int(len(countdict) / 60)
    figure_queue = Queue(6)
    num_img = 1
    sref = sorted(countdict)
    nprocess = 0
    nb_pages = 0
    previous_page = 0
    for pages in range(1,nb_images+2):
        current_page = pages*60
        if nprocess <4:
            nprocess += 1
            process_list.append(Process(target=plot_img,args=(countdict,
                                                              sref[previous_page:current_page],
                                                              figure_queue,)))
            process_list[len(process_list)-1].start()
        else:
            p = process_list.pop(0)
            while figure_queue.empty():
                p.join(1)
            if not figure_queue.empty():
                pdf.savefig(figure_queue.get())
            process_list.append(Process(target=plot_img,args=(countdict,
                                                              sref[previous_page:current_page],
                                                              figure_queue,)))
            process_list[len(process_list)-1].start()
        previous_page = current_page
    while nb_pages < len(process_list):
        if not figure_queue.empty():
            nb_pages += 1
            pdf.savefig(figure_queue.get())
    pdf.close()
    for p in process_list:
        p.join()

def plot_img_par(countdict, ref_list, figure_queue):
    fig, axarr = plt.subplots(15, 4, figsize=(10, 20), sharey='row')
    xax = 0
    yax = 0
    for ref in ref_list:
        m = max(countdict[ref])
        max_coord = len(countdict[ref])
        if m > 0:
            l = [x/m for x in countdict[ref]]
        else:
            l = countdict[ref]
        l_coord = [x/max_coord for x in range(0,max_coord)]
        axarr[xax,yax].plot(l_coord,l)
        axarr[xax,yax].set_title(ref)
        axarr[xax,yax].set_ylim(-0.01,1.1)
        yax += 1
        if yax == 4:
            yax = 0
            xax += 1
    plt.setp([a.get_yticklabels() for a in fig.axes[1:3:60]], visible=True)
    plt.setp([a.set_xticks([]) for a in fig.axes[:-4]])
    plt.setp([a.get_xticklabels() for a in fig.axes[-4:60:2]], visible=False)
    plt.tight_layout(h_pad=0,w_pad=-1.1)
    figure_queue.put(fig)

def plot_coverage(countdict, outfile):
    pdf = PdfPages(outfile)
    pages_number = int(len(countdict)/60)
    uncomplete_last = len(countdict) % 60
    if uncomplete_last > 0:
        pages_number += 1
    ref_list = sorted(countdict)
    line_list = list()
    max_subplots = 60
    for page in range(0, pages_number):
        plot_number = 0
        if page == 0:
            fig, axarr = plt.subplots( 15, 4, figsize=(10,15))
            ax = fig.axes
        for i in range(0, max_subplots):
            try:
                ref = ref_list.pop(0)
            except IndexError:
                if page == pages_number -1:
                    pass
            plot_coverage_plotting(ref, countdict, ax[plot_number],
                                   first = True if page == 0 else False,
                                   line_list=line_list)
            format_figure(fig, ax[plot_number], plot_number,
                          uncomplete_last = uncomplete_last if page == pages_number-1 else 0)
            plot_number += 1
        fig.text(0.02, 0.5, 'Normalized counts', va='center', rotation='vertical', fontsize=12)
        fig.text(0.4, 0.02, 'Normalized coordinates', va='center', rotation='horizontal', fontsize=12)
        fig.suptitle('miRNA coverage maps',fontsize=14)
        fig.tight_layout(rect=[0.05, 0.03, 1, 0.95],h_pad=0.01,w_pad=-1.1)
        pdf.savefig(fig)
    pdf.close()

def format_figure(fig, a, counter, uncomplete_last=0):
    """
    Takes a figure a subplot its position in the figure and the number of
    subplots contained in the last figure.
    """
    if uncomplete_last and counter >= uncomplete_last:
        fig.delaxes(a)
    else:
        if counter % 4 != 0:
            """ If the subplot is not at the left-most position hide y axis """
            plt.setp(a.get_yticklabels(), visible=False)
            plt.setp(a.set_yticks([]))
        elif (counter % 4) % 2 != 0 :
            plt.setp(a.get_yticklabels(), visible=False)
        if counter < 56:
            """ If the subplot isn't at the bottom of the figure hide x axis"""
            plt.setp(a.get_xticklabels(), visible=False)
            plt.setp(a.set_xticks([]))
        elif counter % 2 != 0:
            """ If the subplot is at the bottom of the figure hide x axis for
            every other subplot """
            plt.setp(a.get_xticklabels(), visible=False)

def plot_coverage_plotting(ref, countdict, ax, first=False,
                           last=False, line_list=None,plot_number=None):
    m = max(countdict[ref])
    max_coord = len(countdict[ref])
    if m > 0:
        l = [x/m for x in countdict[ref]]
    else:
        l = countdict[ref]
    l_coord = [x/max_coord for x in range(0,max_coord)]
    ax.set_title(ref, fontsize=7, y=0.91)
    if first:
        line, = ax.plot(l_coord,l)
        line_list.append(line)
        ax.set_ylim(-0.01,1.1)
        return line
    else:
        try:
            line_list[plot_number].set_ydata(l)
            line_list[plot_number].set_xdata(l_coord)
        except TypeError:
            pass

def write_counts(countdict, outfile):
    """
    Takes a dict[<gene_type name>]=count and
    writes a count table
    """
    F = open(outfile, 'w')
    for gene in sorted(countdict):
        F.write('%s\t%s\n' % (gene, str(countdict[gene])))
    F.close()


def main():
    args = Parser()
    bamfile = pysam.AlignmentFile(args.alignment_file, 'rb', check_sq=False)
    if args.pre_mirs:
        pre_mirs = get_pre_mir_counts(bamfile)
        write_counts(pre_mirs, args.pre_mirs)
        if args.pdf_output:
            pre_mirs_coverage = get_pre_mir_coverage(bamfile,
                                                     args.quality_threshold)
            plot_coverage(pre_mirs_coverage, args.pdf_output)
            #write_dataframe_coverage(pre_mirs_coverage, args.lattice)
    if args.mirs:
        mirs = get_mir_counts(bamfile, args.gff_file)
        write_counts(mirs, args.mirs)


if __name__ == '__main__':
    main()
