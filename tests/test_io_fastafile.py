"""Test FastaFile
"""
from basevar.io.fasta import FastaFile


fastafile = "./data/hg19.NC_012920.fasta"
indexfile = "./data/hg19.NC_012920.fasta.fai"


def test_FastaFile(infile):
    fa = FastaFile(infile, infile + ".fai")

    # These are the two inner functions which can just be called by Cython code
    # print("chrM:1 %s" % fa.get_character("chrM", 0))
    # print("chrM:1-10 %s" % fa.get_sequence("chrM", 0, 10))

    print("chr1 length: %d" % fa.get_reference_length("chr1"))
    print("chr2 length: %d" % fa.get_reference_length("chr2"))
    print("chrX length: %d" % fa.get_reference_length("chrX"))
    print("chrY length: %d" % fa.get_reference_length("chrY"))

    print("The total size of reference is: %d" % fa.get_total_sequence_length())
    print("fa filename: %s" % fa.filename)
    print("Number of references seq : %d" % fa.nreferences)
    print("Tuple of reference seq: ", fa.refnames)
    print("Tuple of reference length: ", fa.lengths)

    for c, s in zip(fa.refnames, fa.lengths):
        print("%s length: %d" % (c, s))


if __name__ == "__main__":
    test_FastaFile(fastafile)
