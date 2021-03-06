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

cdef list load_bamdata(dict bamfiles, list samples, bytes chrom, long int start, long int end,
                       char* refseq, options):
    """
    Take a list of BAM files, and a genomic region, and reuturn a list of buffers, containing the
    reads for each BAM file in that region.
    
    ``bamfiles`` is a dict of the object of bamfile path for ``samples``(one for each)
    
    Format in ``bamfiles``: {sample: filename}, it could be create by one line code: 
    bam_objs = {s: f for s, f in zip(samples, input_bamfiles)}
    
    ``samples`` must be the same size as ``bam_objs`` and ``samples`` is the order of input bamfiles
    
    This function could just work for unique sample with only one BAM file. You should merge your 
    bamfiles first if there are multiple BAM files for one sample.
    """

    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead* the_read

    cdef bint is_compress_read = options.is_compress_read
    cdef int qual_bin_size = options.qual_bin_size
    cdef int max_read_thd = options.max_reads

    cdef int total_reads = 0
    cdef int sample_num = len(samples)
    cdef BamReadBuffer sample_read_buffer

    cdef list population_read_buffers = []

    region = "%s:%s-%s" % (chrom, start, end)
    cdef int i
    for i in range(sample_num):
        # assuming the sample is already unique in ``samples``

        reader = Samfile(bamfiles[samples[i]])
        reader.open("r", True)

        # set initial size for BamReadBuffer
        sample_read_buffer = BamReadBuffer(chrom, start, end, options)
        sample_read_buffer.sample = samples[i]

        try:
            reader_iter = reader.fetch(region)
        except Exception as e:
            logger.warning(e.message)
            logger.warning("No data could be retrieved for sample %s in file %s in "
                           "region %s" % (samples[i], reader.filename, region))

            population_read_buffers.append(sample_read_buffer)
            continue

        while reader_iter.cnext():

            the_read = reader_iter.get(0, NULL)
            # if is_compress_read:
            #     compress_read(the_read, refseq, start, end, qual_bin_size)

            sample_read_buffer.add_read_to_buffer(the_read)

            total_reads += 1
            if total_reads > max_read_thd:
                logger.error("Too many reads (%s) in region %s. Quitting now. Either reduce --buffer-size or "
                             "increase --max_reads." % (total_reads, region))

                reader.close()
                sys.exit(1)

            # Todo: we skip all the broken mate reads here, it's that necessary or we should keep them for assembler?

        reader.close()

        # ``population_read_buffers`` will keep the same order as ``samples``,
        # which means will keep the same order as input.
        population_read_buffers.append(sample_read_buffer)

    # return buffers as the same order of input samples/bamfiles
    return population_read_buffers

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

