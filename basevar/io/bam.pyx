"""BAMfile IO"""
import os
import sys

from cpython cimport bool

from basevar.log import logger
from basevar.io.htslibWrapper cimport Samfile


cdef bool is_indexable(filename):
    return filename.lower().endswith((".bam", ".cram"))


cpdef list get_sample_names(list bamfiles, bool filename_has_samplename):
    """Getting sample name in BAM/CRMA files from RG tag and return."""

    logger.info("Start getting all the name of samples from alignment files.")
    if filename_has_samplename:
        logger.info("getting sample name by filename because you set "
                    "'--filename-has-samplename'")

    cdef int file_num = len(bamfiles)
    cdef list sample_names = []
    cdef bytes filename
    for i, al in enumerate(bamfiles):

        if i % 1000 == 0 and i > 0:
            logger.info("loading %d/%d alignment files ..." % (i+1, file_num))

        if filename_has_samplename:
            filename = os.path.basename(al)

            # sample id should be the first element separate by ".",
            # e.g: "CL100045504_L02_61.sorted.rmdup.realign.BQSR.bam", "CL100045504_L02_61" is sample id.
            sample_names.append(filename.split(".")[0])

        else:

            # This may take a very long time to get sampleID from BAM header.
            if is_indexable(al):
                logger.error("Input file %s is not a BAM/CRAM file" % al)
                raise StandardError, "Input file %s is not a BAM/CRAM file" % al

            bf = Samfile(al)
            bf._open("r", True)

            try:

                # header info need to be calculate in bf.header function.
                the_header = bf.header
                if len(the_header["RG"]) > 1:
                    logger.debug("Found multiple read group tags in file %s" % al)

                if "RG" not in the_header["RG"]:
                    logger.error("%s: missing @RG in the header." % al)
                    bf.close()
                    raise StandardError, ("%s: missing @RG in the header." % al)

                sample_names.append(the_header['RG'][0]['SM'])
                del the_header

            except StandardError, e:
                logger.error("Error in BAM header sample parsing. Error was \n%s\n" % e)
                bf.close()
                sys.exit(1)

            bf.close()

    logger.info("Finish loading all %d samples' name\n" % file_num)
    return sample_names


cdef list load_bamdata(list bamfiles, bytes chrom, int start, int end, char* refseq):
    """
    Take a list of BAM files, and a genomic region, and reuturn a list of buffers, containing the
    reads for each BAM file in that region.
    """
    cdef list h = []

    return h















