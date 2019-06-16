"""BAMfile IO"""
import os
import sys


from basevar.log import logger
from basevar.io.read cimport BamReadBuffer
from basevar.io.htslibWrapper cimport Samfile, ReadIterator, cAlignedRead
from basevar.io.htslibWrapper cimport compress_read


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
            if not is_indexable(al):
                logger.error("Input file %s is not a BAM/CRAM file" % al)
                raise StandardError, "Input file %s is not a BAM/CRAM file" % al

            bf = Samfile(al) # load header
            bf._open("r", True)  # load_index

            try:

                # header info need to be calculate in bf.header function.
                the_header = bf.header
                if len(the_header["RG"]) > 1:
                    logger.debug("Found multiple read group tags in file %s" % al)

                if "RG" not in the_header:
                    logger.error("%s: missing @RG in the header." % al)
                    bf.close()
                    raise StandardError, ("%s: missing @RG in the header." % al)

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


cdef list load_bamdata(dict bam_objs, list samples, bytes chrom, long long int start, long long int end,
                       char* refseq, options):
    """
    Take a list of BAM files, and a genomic region, and reuturn a list of buffers, containing the
    reads for each BAM file in that region.
    
    ``bam_objs`` is a dict of the object of class ``Samfile`` for ``samples``(one for each), and the handle 
    has been open to set caching for BAM/CRAM to reduce file IO cost. 
    
    Format in ``bam_objs``: {sample: Samfile}, it could be create by one line code: 
    bam_objs = {s: Samfile(f) for s, f in zip(samples, input_bamfiles)}
    
    ``samples`` must be the same size as ``bam_objs`` and ``samples`` is the order of input bamfiles
    
    This function could just work for unique sample with only one BAM file. You should merge your 
    bamfile first if there are multiple BAM files for one sample.
    """
    cdef list population_read_buffers = []

    cdef bytes sample
    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead* the_read
    cdef int is_compress_read = options.is_compress_read
    cdef int qual_bin_size = options.qual_bin_size
    cdef int max_read_thd = options.max_reads
    cdef int total_reads = 0
    cdef BamReadBuffer sample_read_buffer

    region = "%s:%s-%s" % (chrom, start, end)
    # assuming the sample is already unique in ``samples``
    for i, sample in enumerate(samples):
        assert sample in bam_objs, logger.error("Something is screwy here.")
        reader = bam_objs[sample]

        # Need to lock the thread here when sharing BAM files
        if reader.lock is not None:
            reader.lock.acquire()

        # set initial size for BamReadBuffer
        sample_read_buffer = BamReadBuffer(chrom, start, end, options)
        sample_read_buffer.sample = bytes(sample)
        sample_read_buffer.sample_order = i

        try:
            reader_iter = reader.fetch(region)
        except Exception, e:
            logger.warning(e.message)
            logger.warning("No data could be retrieved for sample %s in file %s in "
                           "region %s" % (sample, reader.filename, region))

            population_read_buffers.append(sample_read_buffer)
            continue

        while reader_iter.cnext():

            the_read = reader_iter.get(0, NULL)
            if is_compress_read:
                compress_read(the_read, refseq, start, end, qual_bin_size)

            sample_read_buffer.add_read_to_buffer(the_read)

            total_reads += 1
            if total_reads % 200000 == 0:
                logger.info("Loaded %s reads in regions %s" % (total_reads, region))

            if total_reads > max_read_thd:
                logger.error("Too many reads (%s) in region %s. Quitting now. Either reduce --buffer-size or "
                             "increase --max_reads." % (total_reads, region))
                for f in bam_objs.values():
                    f.close()

                sys.exit(1)

            # Todo: we skip all the broken mate reads here, it's that necessary or we should keep them for assembler?

        reader.clear_header() # release header to save memory
        # ``pop_read_buffers`` will keep the same order as ``samples``,
        # which means will keep the same order as input.
        population_read_buffers.append(sample_read_buffer)
        # Need to release thread lock here when sharing BAM files
        if reader.lock is not None:
            reader.lock.release()

    # Todo: Do we need to set output order by the order of samples' name?
    cdef list sorted_population_buffers = []
    for sample_read_buffer in population_read_buffers:
        if sample_read_buffer.reads.get_size() > 0:
            sample_read_buffer.chrom_id = sample_read_buffer.reads.array[0].chrom_id

        if not sample_read_buffer.is_sorted:
            sample_read_buffer.sort_reads()

        sample_read_buffer.log_filter_summary()
        sorted_population_buffers.append(sample_read_buffer)

    # return buffers as the same order of input bamfiles
    return sorted_population_buffers















