# cython: profile=True
"""BAMfile IO
"""
import os
import sys

from basevar.log import logger
from basevar.utils cimport c_max
from basevar.io.read cimport BamReadBuffer
from basevar.io.htslibWrapper cimport Samfile, ReadIterator, cAlignedRead
from basevar.caller.batch cimport BatchGenerator

cdef bint is_indexable(filename):
    return filename.lower().endswith((".bam", ".cram"))


cdef list get_sample_names(list bamfiles, bint filename_has_samplename):
    """Getting sample name in BAM/CRMA files from RG tag and return."""

    logger.info("Getting all the samples' name.")
    if filename_has_samplename:
        logger.info("getting sample name by filename because you set "
                    "'--filename-has-samplename'")

    cdef int file_num = len(bamfiles)
    cdef list sample_names = []
    cdef bytes filename
    cdef dict the_header
    cdef int i = 0
    for i in range(file_num):

        if i % 1000 == 0 and i > 0:
            logger.info("loading %d/%d alignment files ..." % (i, file_num))

        if filename_has_samplename:
            filename = os.path.basename(bamfiles[i])

            # sample id should be the first element separate by ".",
            # e.g: "CL100045504_L02_61.sorted.rmdup.realign.BQSR.bam", "CL100045504_L02_61" is sample id.
            sample_names.append(filename.split(".")[0])

        else:

            # This may take a very long time to get sampleID from BAM header if there's a lot of bamfile.
            if not is_indexable(bamfiles[i]):
                logger.error("Input file %s is not a BAM or CRAM file" % bamfiles[i])
                sys.exit(1)

            bf = Samfile(bamfiles[i])
            bf.open("r", True)  # load_index

            try:

                # header info need to be calculate in bf.header function.
                the_header = bf.header
                if len(the_header["RG"]) > 1:
                    logger.debug("Found multiple read group tags in file %s" % bamfiles[i])

                if "RG" not in the_header:
                    logger.error("%s: missing @RG in the header." % bamfiles[i])
                    bf.close()
                    raise StandardError, ("%s: missing @RG in the header." % bamfiles[i])

                sample_names.append(the_header['RG'][0]['SM'])
                bf.clear_header()

            except StandardError, e:
                logger.error("Error in BAM header sample parsing. The error is\n%s\n" % e)
                bf.close()
                sys.exit(1)

            bf.close()

    logger.info("Finish loading all %d samples' names\n" % file_num)
    return sample_names

cdef bint load_data_from_bamfile(Samfile bam_reader,
                                 bytes sample_id,
                                 bytes chrom,
                                 long int start, # 1-base
                                 long int end,   # 1-base
                                 BatchGenerator sample_batch_buffers,
                                 int sample_index,
                                 options):
    """
    This function could just work for unique sample with only one BAM file. You should merge your 
    bamfiles first if there are multiple BAM files for one sample.
    """
    # sample_index is a index label
    if sample_index >= sample_batch_buffers.sample_size:
        logger.error("Index overflow! Index (%d) is lager or equal to sample_size(%d)" % (
            sample_index, sample_batch_buffers.sample_size))
        sys.exit(1)

    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead *the_read

    cdef int total_reads = 0
    cdef BamReadBuffer sample_read_buffer

    cdef long int r_start = c_max(0, start-1)  # make 0-base
    cdef long int r_end = c_max(0, end-1)      # make 0-base
    cdef basestring region = "%s:%s-%s" % (chrom, r_start, r_end)

    # set initial size for BamReadBuffer
    sample_read_buffer = BamReadBuffer(chrom, r_start, r_end, options)
    sample_read_buffer.sample = sample_id

    cdef bint is_empty = True
    try:
        reader_iter = bam_reader.fetch(region)
    except Exception as e:
        logger.warning(e.message)
        logger.warning("No data could be retrieved for sample %s in file %s in "
                       "region %s" % (sample_id, bam_reader.filename, region))
        return is_empty

    while reader_iter.cnext():
        # loading data for one sample in target region
        the_read = reader_iter.get(0, NULL)
        sample_read_buffer.add_read_to_buffer(the_read)
        total_reads += 1

    is_empty = False
    # get batch information for each sample in [start, end]
    sample_batch_buffers.create_batch_in_region(
        (chrom, start, end),
        sample_read_buffer.reads.array,  # this start pointer will move automatically
        sample_read_buffer.reads.array + sample_read_buffer.reads.get_size(),
        sample_index  # ``sample_index`` is sample_index of ``BatchGenerator``
    )

    if options.verbosity > 1:
        logger.info("We get %d good reads for %s from %s" % (total_reads, region, bam_reader.filename))

    return is_empty

