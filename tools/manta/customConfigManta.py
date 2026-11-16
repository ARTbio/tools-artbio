"""Write the configManta.py.ini file used by the Manta Galaxy wrapper.

Settings come from the command line (the "Customize" form) or, with
--configFile, from the [manta] section of a user supplied ini file.
Only the numeric settings listed in SETTINGS are written: any other key
or section of the user file is dropped, and every value must parse as an
integer.

Runs under Python 2 (Manta 1.6 is Python 2) as well as Python 3.
"""
import argparse
import sys

try:
    from ConfigParser import Error as ConfigParserError
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import Error as ConfigParserError
    from configparser import RawConfigParser


MANTA_SECTION = 'manta'

# (name, Manta default, help), in the order they are written to the ini file
SETTINGS = [
    ('minCandidateVariantSize', 8,
     "Run Manta reporting for all SVs/indels at or above this size"),
    ('rnaMinCandidateVariantSize', 1000,
     "Separate option (to provide different default) used for runs in "
     "RNA-mode"),
    ('minEdgeObservations', 3,
     "Remove all edges from the graph unless they're supported by this "
     "many 'observations'"),
    ('graphNodeMaxEdgeCount', 10,
     "If both nodes of an edge have an edge count higher than this, then "
     "skip evaluation of the edge"),
    ('minCandidateSpanningCount', 3,
     "Run discovery and candidate reporting for all SVs/indels with at "
     "least this many spanning support observations"),
    ('minScoredVariantSize', 50,
     "After candidate identification, only score and report SVs/indels at "
     "or above this size"),
    ('minDiploidVariantScore', 10,
     "minimum VCF QUAL score for a variant to be included in the diploid "
     "vcf"),
    ('minPassDiploidVariantScore', 20,
     "VCF QUAL score below which a variant is marked as filtered in the "
     "diploid vcf"),
    ('minPassDiploidGTScore', 15,
     "minimum genotype quality score below which single samples are "
     "filtered for a variant in the diploid vcf"),
    ('minSomaticScore', 10,
     "somatic quality scores below this level are not included in the "
     "somatic vcf"),
    ('minPassSomaticScore', 30,
     "somatic quality scores below this level are filtered in the somatic "
     "vcf"),
    ('enableRemoteReadRetrievalForInsertionsInGermlineCallingModes', 1,
     "GermlineCallingModes includes all calling modes other than "
     "CancerCallingModes"),
    ('enableRemoteReadRetrievalForInsertionsInCancerCallingModes', 0,
     "CancerCallingModes includes tumor-normal subtraction and tumor-only "
     "calling"),
    ('useOverlapPairEvidence', 0,
     "Set 1 if an overlapping read pair will be considered as evidence. "
     "Set to 0 to skip overlapping read pairs"),
]
SETTING_NAMES = [name for name, _, _ in SETTINGS]


def Parser():
    the_parser = argparse.ArgumentParser()
    for name, default, help_text in SETTINGS:
        the_parser.add_argument('--%s' % name, type=int, default=default,
                                help=help_text)
    the_parser.add_argument(
        '--configFile',
        help="Read the settings from the [manta] section of this ini file. "
             "Unknown settings and sections are ignored.")
    return the_parser.parse_args()


def read_config_file(path):
    """Return {setting: int} for the known settings found in path."""
    config = RawConfigParser()
    config.optionxform = str  # keep the case of keys, as Manta does
    try:
        with open(path) as handle:
            if hasattr(config, 'read_file'):
                config.read_file(handle)
            else:
                config.readfp(handle)
    except (ConfigParserError, UnicodeDecodeError) as e:
        sys.exit("Could not parse the configuration file: %s" % e)
    if config.defaults():
        sys.stderr.write("Ignoring section [DEFAULT] of the configuration "
                         "file\n")
        # items() would otherwise merge [DEFAULT] keys into [manta]
        config.defaults().clear()
    for section in config.sections():
        if section != MANTA_SECTION:
            sys.stderr.write("Ignoring section [%s] of the configuration "
                             "file\n" % section)
    values = {}
    if not config.has_section(MANTA_SECTION):
        sys.stderr.write("No [%s] section in the configuration file, using "
                         "default settings\n" % MANTA_SECTION)
        return values
    for key, value in config.items(MANTA_SECTION):
        if key not in SETTING_NAMES:
            sys.stderr.write("Ignoring unsupported setting '%s'\n" % key)
            continue
        try:
            values[key] = int(value)
        except ValueError:
            sys.exit("Setting '%s' must be an integer" % key)
    return values


def write_config(values, path='configManta.py.ini'):
    with open(path, 'w') as handle:
        handle.write('[%s]\n' % MANTA_SECTION)
        handle.write('referenceFasta = /dummy/path/to/genome.fa\n')
        for name in SETTING_NAMES:
            handle.write('%s = %d\n' % (name, values[name]))


if __name__ == "__main__":
    args = Parser()
    values = dict((name, getattr(args, name)) for name in SETTING_NAMES)
    if args.configFile:
        values.update(read_config_file(args.configFile))
    write_config(values)
