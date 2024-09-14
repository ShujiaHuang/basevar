# cython: profile=True
"""BAMfile IO
"""
import os

from libc.stdio cimport fprintf, stderr, stdout
from libc.stdlib cimport exit, EXIT_FAILURE, free
from libc.time cimport clock_t, clock, CLOCKS_PER_SEC

from basevar.io.read cimport BamReadBuffer
from basevar.io.htslibWrapper cimport Samfile, ReadIterator, cAlignedRead
from basevar.caller.batch cimport BatchGenerator

from basevar.utils cimport BaseTypeCmdOptions
from basevar.datatype.strarray cimport StringArray, strarray_init, strarray_append
from basevar.datatype.genomeregion cimport GenomeRegion, make_region_str

cdef bint is_indexable(filename):
    return filename.lower().endswith((".bam", ".cram"))

cdef StringArray get_sample_names(StringArray bamfiles, bint filename_has_samplename):
    """Getting sample name in BAM/CRMA files from RG tag and return."""

    cdef clock_t start_time = clock()

    fprintf(stdout, "[INFO] Getting all the samples' name.\n")
    if filename_has_samplename:
        fprintf(stdout, "[INFO] getting sample name by filename because you set "
                        "'--filename-has-samplename'\n")

    cdef StringArray sample_names
    strarray_init(&sample_names, bamfiles.size)

    cdef str filename
    cdef dict the_header
    cdef int i = 0
    for i in range(bamfiles.size):

        if i % 1000 == 0 and i > 0:
            fprintf(stdout, "[INFO] loading %d/%d alignment files ...\n", i, bamfiles.size)

        if filename_has_samplename:
            filename = os.path.basename(str(bamfiles.array[i]))

            # sample id should be the first element separate by ".",
            # e.g: "CL100045504_L02_61.sorted.rmdup.realign.BQSR.bam", "CL100045504_L02_61" is sample id.
            strarray_append(&sample_names, filename.split(".")[0])

        else:

            # This may take a very long time to get sampleID from BAM header if there's a lot of bamfile.
            if not is_indexable(bamfiles.array[i]):
                fprintf(stderr, "[ERROR] Input file %s is not a BAM or CRAM file\n" % bamfiles.array[i])
                exit(EXIT_FAILURE)

            bf = Samfile(bamfiles.array[i])
            bf.open("r", True)  # load_index

            try:

                # header info need to be calculate in bf.header function.
                the_header = bf.header
                if len(the_header["RG"]) > 1:
                    fprintf(stdout, "[DEBUG] Found multiple read group tags in file %s\n", bamfiles.array[i])

                if "RG" not in the_header:
                    fprintf(stderr, "[ERROR] %s: missing @RG in the header.\n", bamfiles.array[i])
                    bf.close()
                    exit(EXIT_FAILURE)

                strarray_append(&sample_names, the_header['RG'][0]['SM'])
                bf.clear_header()

            except StandardError, e:
                fprintf(stderr, "[ERROR] Error in BAM header sample parsing.\n")
                bf.close()
                exit(EXIT_FAILURE)

            bf.close()

    fprintf(stdout, "[INFO] Finish loading all %ld samples' names, %.1f seconds elapsed\n", bamfiles.size,
            <double>(clock() - start_time)/CLOCKS_PER_SEC)

    return sample_names

cdef list load_bamdata(const StringArray *bamfiles,
                       const StringArray *samples,
                       const GenomeRegion region,  # coordinate in `region` has been make to be 0-base system
                       BaseTypeCmdOptions options):
    """
    Take a list of BAM files, and a genomic region, and reuturn a list of buffers, containing the
    reads for each BAM file in that region.
    
    This function could just work for unique sample with only one BAM file. You should merge your 
    bam first if there are multiple BAM files for one sample.
    """
    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead *the_read
    cdef BamReadBuffer sample_read_buffer

    cdef list population_read_buffers = []
    cdef char *region_string = make_region_str(region)

    cdef unsigned i = 0
    cdef unsigned long total_reads = 0
    for i in range(bamfiles.size):
        # assuming the sample is already unique in ``samples``
        reader = Samfile(bamfiles.array[i])
        reader.open("r", True)

        # set initial size for BamReadBuffer
        sample_read_buffer = BamReadBuffer(samples.array[i], region, options)
        try:
            reader_iter = reader.fetch(region_string)
        except Exception as e:
            fprintf(stdout, "No data could be retrieved for sample %s in file %s in "
                            "region %s\n", samples.array[i], reader.filename, region_string)

            population_read_buffers.append(sample_read_buffer)
            continue

        while reader_iter.cnext():

            the_read = reader_iter.get(0, NULL)
            # if options.is_compress_read:
            #     compress_read(the_read, refseq, start, end, options.qual_bin_size)

            sample_read_buffer.add_read_to_buffer(the_read)
            total_reads += 1
            if total_reads > options.max_reads:
                fprintf(stdout,
                        "[INFO] Too many reads (%ld) in region %s. Quitting now. Either reduce --buffer-count or "
                        "increase --max_reads.\n", total_reads, region_string)

                reader.close()
                exit(EXIT_FAILURE)

            # Todo: I skip all the broken mate reads here, it's that necessary or we should keep them for assembler?

        reader.close()

        # ``population_read_buffers`` will keep the same order as ``samples``, which means will
        # keep the same order as input.
        population_read_buffers.append(sample_read_buffer)

    free(region_string)  # release the memory of 'region' which init in function "utils.set_genome_region"
    # return buffers as the same order of input samples/bamfiles
    return population_read_buffers

cdef bint load_data_from_bamfile(Samfile bam_reader,
                                 char *sample_id,
                                 GenomeRegion region, # coordinate in `region` has been make to be 0-base system
                                 BatchGenerator sample_batch_buffers,
                                 unsigned long sample_index,
                                 options):
    """
    This function just work for unique sample with only one BAM file. You should merge your 
    BAM files if there are multiple BAM files for one sample.
    """
    # sample_index is a index label
    if sample_index >= sample_batch_buffers.sample_size:
        fprintf(stderr, "[ERROR] Index overflow! Index (%d) is lager or equal to sample_size(%d)\n",
                sample_index, sample_batch_buffers.sample_size)
        exit(EXIT_FAILURE)

    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead *the_read

    cdef unsigned long total_reads = 0
    cdef BamReadBuffer sample_read_buffer
    cdef char *region_s = make_region_str(region)

    # set initial size for BamReadBuffer
    sample_read_buffer = BamReadBuffer(sample_id, region, options)

    cdef bint is_empty = True
    try:
        reader_iter = bam_reader.fetch(region_s)
    except Exception as e:
        fprintf(stdout, "[WARNING] No data could be retrieved for sample %s in file %s in "
                       "region %s\n", sample_id, bam_reader.filename, region_s)
        return is_empty

    while reader_iter.cnext():
        # loading data for one sample in target region
        the_read = reader_iter.get(0, NULL)
        sample_read_buffer.add_read_to_buffer(the_read)
        total_reads += 1

    is_empty = False
    # get batch information for each sample in [start, end]
    sample_batch_buffers.create_batch_in_region(
        region,
        sample_read_buffer.reads.array,  # this start pointer will move automatically
        sample_read_buffer.reads.array + sample_read_buffer.reads.get_size(),
        sample_index  # ``sample_index`` is sample_index of ``BatchGenerator``
    )

    if options.verbosity > 1:
        fprintf(stdout, "[INFO] We get %ld good reads in %s from %s\n", total_reads, region, bam_reader.filename)

    free(region_s)
    return is_empty
