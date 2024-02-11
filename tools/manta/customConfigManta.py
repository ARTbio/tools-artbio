import argparse


def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument(
        '--minCandidateVariantSize', type=int, default=8,
        help="Run Manta reporting for all SVs/indels at or above this size")
    the_parser.add_argument(
        '--rnaMinCandidateVariantSize', type=int, default=1000,
        help="Separate option (to provide different default) used for \
              runs in RNA-mode")
    the_parser.add_argument(
        '--minEdgeObservations', type=int, default=3,
        help="Remove all edges from the graph unless they're supported \
              by this many 'observations'")
    the_parser.add_argument(
        '--graphNodeMaxEdgeCount', type=int, default=10,
        help="If both nodes of an edge have an edge count higher than this, \
              then skip evaluation of the edge")
    the_parser.add_argument(
        '--minCandidateSpanningCount', type=int, default=3,
        help="Run discovery and candidate reporting for all SVs/indels with \
              at least this many spanning support observations")
    the_parser.add_argument(
        '--minScoredVariantSize', type=int, default=50,
        help="After candidate identification, only score and report \
              SVs/indels at or above this size")
    the_parser.add_argument(
        '--minDiploidVariantScore', type=int, default=10,
        help="minimum VCF QUAL score for a variant to be included in \
              the diploid vcf")
    the_parser.add_argument(
        '--minPassDiploidVariantScore', type=int, default=20,
        help="VCF QUAL score below which a variant is marked as \
              filtered in the diploid vcf")
    the_parser.add_argument(
        '--minPassDiploidGTScore', type=int, default=15,
        help="minimum genotype quality score below which single samples \
              are filtered for a variant in the diploid vcf")
    the_parser.add_argument(
        '--minSomaticScore', type=int, default=10,
        help="minimum VCF QUAL score for a variant to be included in the \
              diploid vcf")
    the_parser.add_argument(
        '--minPassSomaticScore', type=int, default=30,
        help="somatic quality scores below this level are filtered in the \
              somatic vcf")
    the_parser.add_argument(
        '--enableRemoteReadRetrievalForInsertionsInGermlineCallingModes',
        type=int, default=1,
        help="includes tumor-normal subtraction and tumor-only calling")
    the_parser.add_argument(
        '--enableRemoteReadRetrievalForInsertionsInCancerCallingModes',
        type=int, default=0,
        help="GermlineCallingModes includes all other calling modes")
    the_parser.add_argument(
        '--useOverlapPairEvidence', type=int, default=0,
        help="Set 1 if an overlapping read pair will be considered as \
              evidence. Set to 0 to skip overlapping read pairs")
    args = the_parser.parse_args()
    return args


if __name__ == "__main__":
    args = Parser()
    # recover arguments as a dictionary with keys = argument name and values
    # are argument values
    argsDict = args.__dict__
    ini_lines = []
    # implement first, hard-coded ini lines
    ini_lines.append('[manta]')
    ini_lines.append('referenceFasta = /dummy/path/to/genome.fa')
    # implement the rest of the ini lines for the argsDict
    for argument in argsDict:
        ini_lines.append("%s = %s" % (argument, str(argsDict[argument])))
    # print ini_lines in configManta.py.ini
    handler = open('configManta.py.ini', 'w')
    for line in ini_lines:
        handler.write("%s\n" % line)
