#!/usr/bin/env python
from Bio import SeqIO
import argparse
import sys


def generate_windows(seq, window, step):
    """
    Generates windows of a sequence, with the distance of windows
    defined by *step*.

    seq -- string to split into windows.
    window -- integer specifying the size the generated fragments.
    step -- integer specifiying the distance between adjacent fragments.
    """
    stop = window
    end = len(seq)
    for i in range(stop, end, step):
        start = stop-window
        fragment = seq[start:stop]
        stop_coordinate = stop  # to return real stop coordinate
        stop = stop+step
        yield (fragment, start+1, stop_coordinate)  # start+1 to adjust 0-based range


def write_fragment(description, output_handle, fragment, start, stop):
    """Write out fragments as fasta with description and start/stop coordinates as fasta header"""
    output_string = ">{0}_start:{1}_stop:{2}\n{3}\n".format(description, start, stop, fragment)
    output_handle.write(output_string)


def handle_io(input, output, window=21, step=21):
    """
    Keyword arguments:
    input -- file handle for fasta file containing sequences for which you wish to generate fragments.
    output -- file handle for the multi-fasta that will contain the generated fragments.
    window -- integer specifying the size of the fragments.
    step -- integer specifiying the distance between adjacent fragments.
    """
    record_iterator = SeqIO.parse(input, "fasta")
    for entry in record_iterator:
        seq = str(entry.seq)
        description = str(entry.description)
        windows = generate_windows(seq, window, step)
        [write_fragment(description, output, *fragment) for fragment in windows]
    output.close()
    input.close()


def positive_int(val):
    try:
        assert(int(val) > 0)
    except:
        raise argparse.ArgumentTypeError("'%s' is not a valid positive int" % val)
    return int(val)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fixed size windows in fasta format from multi-fasta sequence.')
    parser.add_argument('--input', type=argparse.FileType('r'), required=True,
                        help='supply an input multi-fasta file.')
    parser.add_argument('--output', type=argparse.FileType('w'), default=sys.stdout,
                        help='supply an output multi-fasta file. If not specified use stdout.')
    parser.add_argument('--window', type=positive_int, default=21,
                        help='Set the size of the generated windows')
    parser.add_argument('--step', type=positive_int, default=21,
                        help='Set distance between the windows')
    args = parser.parse_args()

    handle_io(args.input, args.output, args.window, args.step)
