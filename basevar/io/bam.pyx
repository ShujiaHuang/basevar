"""BAMfile IO"""
import os
import sys

from basevar.log import logger
from basevar.io.read cimport BamReadBuffer
from basevar.io.htslibWrapper cimport Samfile, ReadIterator, cAlignedRead
from basevar.io.htslibWrapper cimport compress_read

cdef int LOW_QUAL_BASES = 0
cdef int UNMAPPED_READ = 1
cdef int MATE_UNMAPPED = 2
cdef int MATE_DISTANT = 3
cdef int SMALL_INSERT = 4
cdef int DUPLICATE = 5
cdef int LOW_MAP_QUAL = 6


cdef extern from "stdlib.h":
    void *malloc(size_t)
    void *calloc(size_t, size_t)
    void free(void *)


cdef bint is_indexable(filename):
    return filename.lower().endswith((".bam", ".cram"))


cdef list get_sample_names(list bamfiles, bint filename_has_samplename):
    """Getting sample name in BAM/CRMA files from RG tag and return."""

    logger.info("Start getting all the name of samples from alignment files.")
    if filename_has_samplename:
        logger.info("getting sample name by filename because you set "
                    "'--filename-has-samplename'")

    cdef int file_num = len(bamfiles)
    cdef list sample_names = []
    cdef bytes filename
    cdef int i = 0
    for i in range(file_num):

        if i % 1000 == 0 and i > 0:
            logger.info("loading %d/%d alignment files ..." % (i+1, file_num))

        if filename_has_samplename:
            filename = os.path.basename(bamfiles[i])

            # sample id should be the first element separate by ".",
            # e.g: "CL100045504_L02_61.sorted.rmdup.realign.BQSR.bam", "CL100045504_L02_61" is sample id.
            sample_names.append(filename.split(".")[0])

        else:

            # This may take a very long time to get sampleID from BAM header.
            if not is_indexable(bamfiles[i]):
                logger.error("Input file %s is not a BAM/CRAM file" % bamfiles[i])
                raise StandardError, "Input file %s is not a BAM/CRAM file" % bamfiles[i]

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
                del the_header

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

        # Need to lock the thread here when sharing BAM files
        if reader.lock is not None:
            reader.lock.acquire()

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
            if is_compress_read:
                compress_read(the_read, refseq, start, end, qual_bin_size)

            sample_read_buffer.add_read_to_buffer(the_read)

            total_reads += 1
            if total_reads > max_read_thd:
                logger.error("Too many reads (%s) in region %s. Quitting now. Either reduce --buffer-size or "
                             "increase --max_reads." % (total_reads, region))

                reader.close()
                sys.exit(1)

            # Todo: we skip all the broken mate reads here, it's that necessary or we should keep them for assembler?

        # close bamfile
        reader.close()

        # ``population_read_buffers`` will keep the same order as ``samples``,
        # which means will keep the same order as input.
        population_read_buffers.append(sample_read_buffer)

        # Need to release thread lock here when sharing BAM files
        if reader.lock is not None:
            reader.lock.release()

    # return buffers as the same order of input samples/bamfiles
    return population_read_buffers

    # cdef list sorted_population_read_buffers = []
    # for sample_read_buffer in population_read_buffers:
    #     if sample_read_buffer.reads.get_size() > 0:
    #         sample_read_buffer.chrom_id = sample_read_buffer.reads.array[0].chrom_id
    #
    #     if not sample_read_buffer.is_sorted:
    #         sample_read_buffer.sort_reads()
    #
    #     log_filter_summary(sample_read_buffer, options.verbosity)
    #     sorted_population_read_buffers.append(sample_read_buffer)
    #
    # # return buffers as the same order of input samples/bamfiles
    # return sorted_population_read_buffers


cdef void log_filter_summary(BamReadBuffer read_buffer, int verbosity):
        """Useful debug information about which reads have been filtered out.
        """
        if verbosity >= 3:
            region = "%s:%s-%s" % (read_buffer.chrom, read_buffer.start+1, read_buffer.end+1)
            logger.debug("Sample %s has %s good reads in %s" % (read_buffer.sample, read_buffer.reads.get_size(), region))
            logger.debug("Sample %s has %s bad reads in %s" % (read_buffer.sample, read_buffer.bad_reads.get_size(), region))
            logger.debug("Sample %s has %s broken mates %s" % (read_buffer.sample, read_buffer.broken_mates.get_size(), region))
            logger.debug("|-- low map quality reads = %s" % (read_buffer.filtered_read_counts_by_type[LOW_MAP_QUAL]))
            logger.debug("|-- low qual reads = %s" % (read_buffer.filtered_read_counts_by_type[LOW_QUAL_BASES]))
            logger.debug("|-- un-mapped reads = %s" % (read_buffer.filtered_read_counts_by_type[UNMAPPED_READ]))
            logger.debug("|-- reads with unmapped mates = %s" % (read_buffer.filtered_read_counts_by_type[MATE_UNMAPPED]))
            logger.debug("|-- reads with distant mates = %s" % (read_buffer.filtered_read_counts_by_type[MATE_DISTANT]))
            logger.debug("|-- reads pairs with small inserts = %s" % (read_buffer.filtered_read_counts_by_type[SMALL_INSERT]))
            logger.debug("|__ duplicate reads = %s\n" % (read_buffer.filtered_read_counts_by_type[DUPLICATE]))

            if read_buffer.trim_overlapping == 1:
                logger.debug("Overlapping segments of read pairs were clipped")
            else:
                logger.debug("Overlapping segments of read pairs were not clipped")

            logger.debug("Overhanging bits of reads were not clipped")













