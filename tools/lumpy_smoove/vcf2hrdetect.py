#!/usr/bin/env python
import argparse
import re
import sys


def create_arg_parser():
    """Creates and returns the argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Convert a VCF file from lumpy-smoove to a tabular format "
            "compatible with the HRDetect pipeline."
        )
    )
    parser.add_argument(
        'vcf_file',
        help='Path to the input VCF file.'
    )
    parser.add_argument(
        'output_file',
        help='Path to the output tabular file.'
    )
    return parser


def parse_breakend_alt(alt_field):
    """
    Parses the ALT field for a breakend and returns chromosome and position.

    Args:
        alt_field (str): The ALT field (column 5) of a VCF line.

    Returns:
        tuple: A tuple containing (chromosome, position) or (None, None)
               if parsing fails.
    """
    # Search for patterns ]chr:pos] or [chr:pos[
    pattern = (
        r"\](?P<chrom1>[^:]+):(?P<pos1>\d+)\]|"
        r"\[(?P<chrom2>[^:]+):(?P<pos2>\d+)\["
    )
    match = re.search(pattern, alt_field)

    if not match:
        return None, None

    groups = match.groupdict()
    chrom = groups['chrom1'] or groups['chrom2']
    pos = groups['pos1'] or groups['pos2']
    return chrom, pos


def process_vcf(vcf_path, output_path):
    """
    Reads a VCF file, converts it, and writes the result to a tabular file.

    Args:
        vcf_path (str): Path to the input VCF file.
        output_path (str): Path to the output tabular file.
    """
    header = ["chr1", "pos1", "chr2", "pos2", "type"]
    try:
        with open(vcf_path, 'r') as infile, open(output_path, 'w') as outfile:
            outfile.write("\t".join(header) + "\n")

            for line in infile:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                chrom1 = fields[0]
                pos1 = fields[1]
                info = fields[7]

                # Attempt to extract the structural variant type from the info
                svtype_match = re.search(r'SVTYPE=([^;]+)', info)
                if not svtype_match:
                    continue  # Skip lines without SVTYPE tag
                svtype = svtype_match.group(1)

                if svtype == "BND":  # Breakend (INV or TRA)
                    alt_field = fields[4]
                    chrom2, pos2 = parse_breakend_alt(alt_field)
                    if not (chrom2 and pos2):
                        continue
                    event_type = "INV" if chrom1 == chrom2 else "TRA"
                    row = [chrom1, pos1, chrom2, pos2, event_type]
                    outfile.write("\t".join(row) + "\n")

                else:  # Other SV types (DEL, DUP, etc.)
                    end_match = re.search(r'END=([^;]+)', info)
                    if not end_match:
                        continue
                    pos2 = end_match.group(1)
                    chrom2 = chrom1
                    row = [chrom1, pos1, chrom2, pos2, svtype]
                    outfile.write("\t".join(row) + "\n")

    except FileNotFoundError:
        print(f"Error: File '{vcf_path}' not found.",
              file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"IO Error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main function of the script."""
    parser = create_arg_parser()
    args = parser.parse_args()
    process_vcf(args.vcf_file, args.output_file)


if __name__ == '__main__':
    main()
