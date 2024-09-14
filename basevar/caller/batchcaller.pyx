# cython: profile=True
"""
Package for parsing bamfile
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
from libc.stdio cimport fprintf, stderr, stdout
from libc.stdlib cimport exit, EXIT_FAILURE, free

from basevar.io.openfile import Open
from basevar.io.fasta cimport FastaFile
from basevar.io.bam cimport load_bamdata
from basevar.io.read cimport BamReadBuffer
from basevar.caller.batch cimport BatchGenerator, BatchInfo

from basevar.utils cimport BaseTypeCmdOptions, c_max
from basevar.datatype.strarray cimport StringArray, convert_strarray_to_string
from basevar.datatype.genomeregion cimport GenomeRegion
from basevar.datatype.dynamicstring cimport dstring, dstring_destroy


cdef void generate_batchfile(const GenomeRegion region,  # 1-base in GenomeRegion
                             const StringArray *batch_align_files,
                             const StringArray *batch_sample_ids,
                             FastaFile fa,
                             char *out_batch_file,
                             BaseTypeCmdOptions options):

    """Loading bamfile and create a batchfile in ``region``.

    Parameters:
        ``region``: The coordinate in regions is 1-base system.
    """
    cdef GenomeRegion region_0_base  # A genome region in 0-base system
    cdef list sample_read_buffers
    try:
        # load the whole mapping reads in [chrom_name, start, end]
        region_0_base.chrom = region.chrom
        region_0_base.start = c_max(region.start-1, 0)
        region_0_base.end = c_max(0, region.end-1)
        sample_read_buffers = load_bamdata(batch_align_files, batch_sample_ids, region_0_base, options)

    except Exception as e:
        fprintf(stderr, "[ERROR] Exception error in region %s:%ld-%ld.\n", region.chrom, region.start, region.end)
        exit(EXIT_FAILURE)

    if sample_read_buffers is None or len(sample_read_buffers) == 0:
        fprintf(stdout, "Skipping region %s:%ld-%ld as it's empty.\n", region.chrom, region.start, region.end)
        return

    cdef BatchGenerator batch_buffer
    batch_buffer = BatchGenerator(region, fa, batch_sample_ids.size, options)

    cdef unsigned sample_index = 0
    cdef unsigned longest_read_size = 0
    cdef BamReadBuffer sample_read_buffer
    # loop all samples
    for sample_index in range(batch_sample_ids.size):
        sample_read_buffer = sample_read_buffers[sample_index]
        if longest_read_size < sample_read_buffer.reads.get_length_of_longest_read():
            longest_read_size = sample_read_buffer.reads.get_length_of_longest_read()

        # get a batch information for each sample in [start, end], start and end are the two element of ``region``
        batch_buffer.create_batch_in_region(
            region,
            sample_read_buffer.reads.array,  # this start pointer will move automatically
            sample_read_buffer.reads.array + sample_read_buffer.reads.get_size(),
            sample_index  # sample_index is the index in ``BatchGenerator``
        )

    # Todo: take care, although this code may not been called forever.
    if longest_read_size > options.r_len:
        options.r_len = longest_read_size

    output_batch_file(batch_buffer, out_batch_file, batch_sample_ids)
    return

cdef void output_batch_file(BatchGenerator batch_buffer, char *out_batch_file, const StringArray *batch_sample_ids):

    cdef int position_number = 0
    cdef int j = 0

    cdef BatchInfo batch_info
    cdef dstring ds

    cdef char *sample_string = convert_strarray_to_string(batch_sample_ids, ",")
    with Open(out_batch_file, "wb", isbgz=True) as OUT:

        OUT.write("##fileformat=BaseVarBatchFile_v1.0\n")
        if batch_sample_ids.size:
            OUT.write("##SampleIDs=%s\n" % sample_string)

        suff_header = ["#CHROM", "POS", "REF", "Depth(CoveredSample)", "MappingQuality", "Readbases",
                       "ReadbasesQuality", "ReadPositionRank", "Strand"]
        OUT.write("%s\n" % "\t".join(suff_header))

        position_number = len(batch_buffer.batch_heap)
        for j in range(position_number):
            batch_info = batch_buffer.batch_heap[j]
            ds = batch_info.get_str()

            OUT.write("%s\n" % str(ds.s))
            dstring_destroy(&ds)

    if sample_string != NULL:
        free(sample_string)

    return
