"""Test bam io
"""
from basevar.io.bam import get_sample_names

base_dir = "./data/140k_thalassemia_brca_bam"
bamfile = "./data/140k_thalassemia_brca_bam/bam90.list"


def test_get_sample_names(bamfiles):

    bfs = get_sample_names(bamfiles, False)
    print(bfs)

    return


if __name__ == "__main__":

    print("Startting ......")

    bamfiles = []
    with open(bamfile) as F:
        for line in F:
            bamfiles.append(base_dir + "/" + line.strip())

    test_get_sample_names(bamfiles)
