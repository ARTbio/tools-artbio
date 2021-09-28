import argparse


description = ("This script will create a configuration file for samples \
                to be run in Pindel")


parser = argparse.ArgumentParser()
parser.add_argument("--input_file", nargs="*",
                    help="One or more alignment files")
parser.add_argument("--insert_size", nargs="+",
                    help="Expected Insert size")
parser.add_argument("--sample_label", nargs="+",
                    help="Sample label")
parser.add_argument("--output_config_file",
                    help="Output config file")
args = parser.parse_args()

template = "{input_file}\t{insert_size}\t{sample_label}\n"
with open(args.output_config_file, "w") as output:
    for input_file, insert_size, sample_label in zip(args.input_file,
                                                     args.insert_size,
                                                     args.sample_label):
        config_line = template.format(input_file=input_file,
                                      insert_size=insert_size,
                                      sample_label=sample_label)
        output.write(config_line)


