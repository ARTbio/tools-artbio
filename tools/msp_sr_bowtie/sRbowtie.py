#!/usr/bin/env python
# small RNA oriented bowtie wrapper
# version 1.5 17-7-2014: arg parser implementation
# Usage sRbowtie.py <1 input_fasta_file> <2 alignment method> <3 -v mismatches> <4 out_type> <5 buildIndexIfHistory> <6 fasta/bowtie index> <7 bowtie output> <8 ali_fasta> <9 unali_fasta> <10 --num-threads \${GALAXY_SLOTS:-4}>
# current rev: for bowtie __norc, move from --supress 2,6,7,8 to --supress 6,7,8. Future Parser must be updated to take into account this standardisation
# Christophe Antoniewski <drosofff@gmail.com>

import sys
import os
import subprocess
import tempfile
import shutil
import argparse


def Parser():
    the_parser = argparse.ArgumentParser(
        description="bowtie wrapper for small fasta reads")
    the_parser.add_argument(
        '--input', action="store", type=str, help="input file")
    the_parser.add_argument(
        '--input-format', dest="input_format", action="store", type=str, help="fasta or fastq")
    the_parser.add_argument('--method', action="store", type=str,
                            help="RNA, unique, multiple, k_option, n_option, a_option")
    the_parser.add_argument('--v-mismatches', dest="v_mismatches", action="store",
                            type=str, help="number of mismatches allowed for the alignments")
    the_parser.add_argument(
        '--output-format', dest="output_format", action="store", type=str, help="tabular, sam, bam")
    the_parser.add_argument(
        '--output', action="store", type=str, help="output file path")
    the_parser.add_argument(
        '--index-from', dest="index_from", action="store", type=str, help="indexed or history")
    the_parser.add_argument('--index-source', dest="index_source",
                            action="store", type=str, help="file path to the index source")
    the_parser.add_argument(
        '--aligned', action="store", type=str, help="aligned read file path, maybe None")
    the_parser.add_argument('--unaligned', action="store",
                            type=str, help="unaligned read file path, maybe None")
    the_parser.add_argument('--num-threads', dest="num_threads",
                            action="store", type=str, help="number of bowtie threads")
    args = the_parser.parse_args()
    return args


def stop_err(msg):
    sys.stderr.write('%s\n' % msg)
    sys.exit()


def bowtieCommandLiner(alignment_method="RNA", v_mis="1", out_type="tabular",
                       aligned="None", unaligned="None", input_format="fasta", input="path",
                       index="path", output="path", pslots="4"):
    if input_format == "fasta":
        input_format = "-f"
    elif (input_format == "fastq") or (input_format == "fastqsanger"):
        input_format = "-q"
    else:
        raise Exception('input format must be one of fasta or fastq')
    if alignment_method == "RNA":
        x = "-v %s -M 1 --best --strata -p %s --norc --suppress 6,7,8" % (
            v_mis, pslots)
    elif alignment_method == "unique":
        x = "-v %s -m 1 -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif alignment_method == "multiple":
        x = "-v %s -M 1 --best --strata -p %s --suppress 6,7,8" % (
            v_mis, pslots)
    elif alignment_method == "k_option":
        x = "-v %s -k 1 --best -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif alignment_method == "n_option":
        x = "-n %s -M 1 --best -p %s --suppress 6,7,8" % (v_mis, pslots)
    elif alignment_method == "a_option":
        x = "-v %s -a --best -p %s --suppress 6,7,8" % (v_mis, pslots)
    if aligned == "None" and unaligned == "None":
        fasta_command = ""
    elif aligned != "None" and unaligned == "None":
        fasta_command = " --al %s" % aligned
    elif aligned == "None" and unaligned != "None":
        fasta_command = " --un %s" % unaligned
    else:
        fasta_command = " --al %s --un %s" % (aligned, unaligned)
    x = x + fasta_command
    if out_type == "tabular":
        return "bowtie %s %s %s %s > %s" % (x, index, input_format, input, output)
    elif out_type == "sam":
        return "bowtie %s -S %s %s %s > %s" % (x, index, input_format, input, output)
    elif out_type == "bam":
        return "bowtie %s -S %s %s %s |samtools view -bS - > %s" % (
            x, index, input_format, input, output)


def bowtie_squash(fasta):
    # make temp directory for bowtie indexes
    tmp_index_dir = tempfile.mkdtemp()
    ref_file = tempfile.NamedTemporaryFile(dir=tmp_index_dir)
    ref_file_name = ref_file.name
    # by default, delete the temporary file, but ref_file.name is now stored
    # in ref_file_name
    ref_file.close()
    # symlink between the fasta source file and the deleted ref_file name
    os.symlink(fasta, ref_file_name)
    # bowtie command line, which will work after changing dir
    # (cwd=tmp_index_dir)
    cmd1 = 'bowtie-build -f %s %s' % (ref_file_name, ref_file_name)
    try:
        FNULL = open(os.devnull, 'w')
        # a path string for a temp file in tmp_index_dir. Just a string
        tmp = tempfile.NamedTemporaryFile(dir=tmp_index_dir).name
        # creates and open a file handler pointing to the temp file
        tmp_stderr = open(tmp, 'wb')
        # both stderr and stdout of bowtie-build are redirected in  dev/null
        proc = subprocess.Popen(
            args=cmd1, shell=True, cwd=tmp_index_dir, stderr=FNULL, stdout=FNULL)
        returncode = proc.wait()
        tmp_stderr.close()
        FNULL.close()
        sys.stdout.write(cmd1 + "\n")
    except Exception as e:
        # clean up temp dir
        if os.path.exists(tmp_index_dir):
            shutil.rmtree(tmp_index_dir)
            stop_err('Error indexing reference sequence\n' + str(e))
    # no Cleaning if no Exception, tmp_index_dir has to be cleaned after
    # bowtie_alignment()
    # bowtie fashion path without extention
    index_full_path = os.path.join(tmp_index_dir, ref_file_name)
    return tmp_index_dir, index_full_path


def bowtie_alignment(command_line, flyPreIndexed=''):
    # make temp directory just for stderr
    tmp_index_dir = tempfile.mkdtemp()
    tmp = tempfile.NamedTemporaryFile(dir=tmp_index_dir).name
    tmp_stderr = open(tmp, 'wb')
    # conditional statement for sorted bam generation viewable in Trackster
    if "samtools" in command_line:
        # recover the final output file name
        target_file = command_line.split()[-1]
        path_to_unsortedBam = os.path.join(tmp_index_dir, "unsorted.bam")
        path_to_sortedBam = os.path.join(tmp_index_dir, "unsorted.bam.sorted")
        first_command_line = " ".join(
            command_line.split()[:-3]) + " -o " + path_to_unsortedBam + " - "
        # example: bowtie -v 0 -M 1 --best --strata -p 12 --suppress 6,7,8 -S
        # /home/galaxy/galaxy-dist/bowtie/Dmel/dmel-all-chromosome-r5.49 -f
        # /home/galaxy/galaxy-dist/database/files/003/dataset_3460.dat
        # |samtools view -bS -o /tmp/tmp_PgMT0/unsorted.bam -
        # generates an "unsorted.bam.sorted.bam file", NOT an
        # "unsorted.bam.sorted" file
        second_command_line = "samtools sort  %s %s" % (
            path_to_unsortedBam, path_to_sortedBam)
        # fileno() method return the file descriptor number of tmp_stderr
        p = subprocess.Popen(
            args=first_command_line, cwd=tmp_index_dir, shell=True, stderr=tmp_stderr.fileno())
        returncode = p.wait()
        sys.stdout.write("%s\n" % first_command_line + str(returncode))
        p = subprocess.Popen(
            args=second_command_line, cwd=tmp_index_dir, shell=True, stderr=tmp_stderr.fileno())
        returncode = p.wait()
        sys.stdout.write("\n%s\n" % second_command_line + str(returncode))
        if os.path.isfile(path_to_sortedBam + ".bam"):
            shutil.copy2(path_to_sortedBam + ".bam", target_file)
    else:
        p = subprocess.Popen(
            args=command_line, shell=True, stderr=tmp_stderr.fileno())
        returncode = p.wait()
        sys.stdout.write(command_line + "\n")
    tmp_stderr.close()
    # cleaning if the index was created in the fly
    if os.path.exists(flyPreIndexed):
        shutil.rmtree(flyPreIndexed)
    # cleaning tmp files and directories
    if os.path.exists(tmp_index_dir):
        shutil.rmtree(tmp_index_dir)
    return


def __main__():
    args = Parser()
    F = open(args.output, "w")
    if args.index_from == "history":
        tmp_dir, index_path = bowtie_squash(args.index_source)
    else:
        tmp_dir, index_path = "dummy/dymmy", args.index_source
    command_line = bowtieCommandLiner(args.method, args.v_mismatches, args.output_format,
                                      args.aligned, args.unaligned, args.input_format, args.input, 
                                      index_path, args.output, args.num_threads)
    bowtie_alignment(command_line, flyPreIndexed=tmp_dir)
    F.close()
if __name__ == "__main__":
    __main__()
