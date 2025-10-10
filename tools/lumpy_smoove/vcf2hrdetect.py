import sys
import argparse
import re


def parse_args():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert a lumpy-smoove VCF to the HRDetect tabular format."
    )
    parser.add_argument(
        "input_vcf", help="Path to the input VCF file from lumpy-smoove."
    )
    parser.add_argument(
        "output_tabular", help="Path for the output tabular file."
    )
    return parser.parse_args()


def parse_info_tag(info_str, tag_name):
    """
    Safely parses a tag from the VCF INFO field.
    Example: SVTYPE=DEL;END=12345 -> parse_info_tag(info, "END") -> "12345"
    """
    # Using regex to find the tag and its value
    match = re.search(f'{tag_name}=([^;]+)', info_str)
    if match:
        return match.group(1)
    return None


def parse_breakend_alt(alt_field):
    """
    Parses the ALT field for a breakend (BND) variant to find the mate's
    coordinates.
    Example: ]chr1:12345]N or N[chr1:12345[
    """
    # Break long regex pattern across multiple lines for readability
    pattern = (
        r'\](?P<chrom1>[^:]+):(?P<pos1>\d+)\]'  # Case 1: ]chr:pos]
        r'|'                                   # OR
        r'\[(?P<chrom2>[^:]+):(?P<pos2>\d+)\[') # Case 2: [chr:pos[
    match = re.search(pattern, alt_field)
    if match:
        # The regex captures into named groups, only one of which will be populated
        groups = match.groupdict()
        chrom = groups['chrom1'] or groups['chrom2']
        pos = groups['pos1'] or groups['pos2']
        return chrom, pos
    return None, None


def process_vcf(input_vcf_path, output_tabular_path):
    """
    Reads the VCF file line by line, processes SV records, and writes to the
    output file.
    """
    try:
        with open(input_vcf_path, 'r') as infile:
            with open(output_tabular_path, 'w') as outfile:
                # Write header
                outfile.write("chr1\tpos1\tchr2\tpos2\ttype\n")

                for line in infile:
                    if line.startswith("#"):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue  # Skip malformed lines

                    chrom1 = fields[0]
                    pos1 = fields[1]
                    sv_id = fields[2]
                    alt = fields[4]
                    info = fields[7]

                    # Case 1: Paired breakend events (Translocations / Inversions)
                    # These are identified by a "_1" suffix in their ID.
                    # The mate's location is in the ALT field.
                    if '_1' in sv_id:
                        chrom2, pos2 = parse_breakend_alt(alt)
                        if not chrom2 or not pos2:
                            continue  # Skip if ALT is not a valid breakend

                        sv_type = "INV" if chrom1 == chrom2 else "TRA"
                        outfile.write(
                            f"{chrom1}\t{pos1}\t{chrom2}\t{pos2}\t{sv_type}\n"
                        )

                    # Case 2: Simple intra-chromosomal events (Dels, Dups, etc.)
                    # These do not have "_" in their ID. End pos is in INFO.
                    elif '_' not in sv_id:
                        sv_type = parse_info_tag(info, "SVTYPE")
                        end_pos = parse_info_tag(info, "END")
                        if not sv_type or not end_pos:
                            continue  # Skip if required INFO tags are missing

                        chrom2 = chrom1  # By definition, intra-chromosomal
                        outfile.write(
                            f"{chrom1}\t{pos1}\t{chrom2}\t{end_pos}\t{sv_type}\n"
                        )

                    # Note: Records with "_2" in ID are implicitly ignored.

    except FileNotFoundError:
        print(
            f"Error: Input file not found at {input_vcf_path}",
            file=sys.stderr
        )
        sys.exit(1)


def main():
    """
    Main function to run the script.
    """
    args = parse_args()
    process_vcf(args.input_vcf, args.output_tabular)


if __name__ == "__main__":
    main()

