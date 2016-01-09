metavisitor

version 1_1:
- introduce the tool blast_to_contig which accomodate both blastn and tblastx alignments
to reconstruct a final genome scaffold with de novo assembled sequences (capitals) and
sequences from the Guide (small letters)
- use the changeset_revision 7f3c448e119b of the ncbi_blast_plus suite
- correct minors bug in names of the files returned by the fetch_fasta_from_ncbi tool