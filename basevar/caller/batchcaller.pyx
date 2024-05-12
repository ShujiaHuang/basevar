# cython: profile=True
"""
Package for parsing bamfile
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
import os
import sys
import time

from basevar.log import logger

from basevar.io.openfile import Open
from basevar.io.fasta cimport FastaFile
from basevar.io.bam cimport load_bamdata
from basevar.io.read cimport BamReadBuffer
from basevar.caller.batch cimport BatchGenerator


cdef list create_batchfiles_in_regions(bytes chrom_name,
                                       list regions,
                                       long int region_boundary_start, # 1-base
                                       long int region_boundary_end,   # 1-base
                                       list align_files,
                                       FastaFile fa,
                                       list samples,
                                       basestring outdir,
                                       object options):
    """
    ``regions`` is a 2-D array
        They all are the some chromosome: [[start1,end1], [start2, end2], ...]
        
    ``samples``: The sample id of align_files
    ``fa``:
        # get sequence of chrom_name from reference fasta
        fa = self.ref_file_hd.fetch(chrid)

    """
    # store all the batch files
    cdef list batchfiles = []
    cdef int batchcount = options.batch_count
    cdef int part_num = int(len(align_files) / batchcount)
    if part_num * batchcount < len(align_files):
        part_num += 1

    cdef bytes refseq_bytes = fa.get_sequence(chrom_name, region_boundary_start,
                                              region_boundary_end+1.0*options.r_len)
    cdef char* refseq = refseq_bytes

    cdef int m = 0, i = 0
    for i in range(0, len(align_files), batchcount):
        # Create a batch of temp files which we call them batchfiles for variant discovery
        start_time = time.time()

        m += 1
        part_file_name = "basevar.%s.%d_%d.batch.gz" % (
            "%s.%s.%s" % (chrom_name, region_boundary_start+1, region_boundary_end+1), m, part_num)
        part_file_name = os.path.join(outdir, part_file_name)  # Join Path could fix different OS

        # store the name of batchfiles into a list.
        batchfiles.append(part_file_name)
        if options.smartrerun and os.path.isfile(part_file_name):
            # ``part_file_name`` is exists We don't have to create it again if setting `smartrerun`
            logger.info("%s already exists, we don't have to create it again, "
                        "when you set `smartrerun`" % part_file_name)
            continue
        else:
            logger.info("Creating batchfile %s\n" % part_file_name)

        # One batch of alignment files, the size and order are the same with ``sub_align_files`` and
        # ``batch_sample_ids``
        sub_align_files = align_files[i:i+batchcount]
        batch_sample_ids = None
        if samples:
            batch_sample_ids = samples[i:i+batchcount]

        generate_batchfile(chrom_name,
                           region_boundary_start,  # 1-base
                           region_boundary_end,    # 1-base
                           regions,
                           sub_align_files,
                           refseq,
                           fa,
                           part_file_name,
                           batch_sample_ids,
                           options)

        logger.info("Done for batchfile %s , %d seconds elapsed." % (
            part_file_name, time.time() - start_time))

    return batchfiles


cdef void generate_batchfile(bytes chrom_name,
                             long int bigstart, # 1-base
                             long int bigend,   # 1-base
                             list regions,
                             list batch_align_files,
                             char *ref_seq,
                             FastaFile fa,
                             bytes out_batch_file,
                             list batch_sample_ids,
                             object options):

    """Loading bamfile and create a batchfile in ``regions``.

    Parameters:
        ``bigstart``: It's already 0-base position
        ``bigend``: It's already 0-base position
        ``regions``: The coordinate in regions is 1-base system.
    """
    # Make bamfiles dict: {sample_id:bamfile,...,}
    cdef int sample_size = len(batch_sample_ids)
    cdef dict bam_files = {}
    cdef int i = 0
    for i in range(sample_size):
        bam_files[batch_sample_ids[i]] = batch_align_files[i]

    cdef list sample_read_buffers
    try:
        # load the whole mapping reads in [chrom_name, bigstart, bigend]
        sample_read_buffers = load_bamdata(bam_files, batch_sample_ids, chrom_name,
                                           bigstart-1, bigend-1, ref_seq, options)

    except Exception, e:
        logger.error("Exception in region %s:%s-%s. Error: %s" % (chrom_name, bigstart, bigend, e))
        sys.exit(1)

    if sample_read_buffers is None or len(sample_read_buffers) == 0:
        logger.info("Skipping region %s:%s-%s as it's empty." % (chrom_name, bigstart, bigend))
        return

    cdef BatchGenerator batch_buffer
    cdef list region_batch_buffers = []

    cdef int sample_index
    cdef int longest_read_size = 0
    cdef BamReadBuffer sample_read_buffer
    cdef long int reg_start, reg_end
    for reg_start, reg_end in regions:
        # initialization the BatchGenerator in `ref_name:reg_start-reg_end`
        batch_buffer = BatchGenerator(chrom_name, reg_start, reg_end, fa, sample_size, options)

        # loop all samples
        for sample_index in range(sample_size):
            sample_read_buffer = sample_read_buffers[sample_index]
            if longest_read_size < sample_read_buffer.reads.get_length_of_longest_read():
                longest_read_size = sample_read_buffer.reads.get_length_of_longest_read()

            # get batch information for each sample in [start, end]
            batch_buffer.create_batch_in_region(
                (chrom_name, reg_start, reg_end),
                sample_read_buffer.reads.array,  # this start pointer will move automatically
                sample_read_buffer.reads.array + sample_read_buffer.reads.get_size(),
                sample_index  # sample_index is the index in ``BatchGenerator``
            )

        # store information in each [reg_start, reg_end] of all samples
        region_batch_buffers.append(batch_buffer)

    # Todo: take care, although this code may not been called forever.
    if longest_read_size > options.r_len:
        options.r_len = longest_read_size

    output_batch_file(chrom_name, fa, region_batch_buffers, out_batch_file, batch_sample_ids, regions)
    return


cdef void output_batch_file(bytes chrom_name, FastaFile fa, list region_batch_buffers, bytes out_batch_file,
                            list batch_sample_ids, list regions):

    cdef int region_number = len(region_batch_buffers)
    cdef int position_number = 0
    cdef int i = 0, j = 0
    cdef BatchGenerator batch_buffer
    with Open(out_batch_file, "wb", isbgz=True) if out_batch_file.endswith(".gz") else \
            open(out_batch_file, "w") as OUT:

        OUT.write("##fileformat=BaseVarBatchFile_v1.0\n")
        if batch_sample_ids:
            OUT.write("##SampleIDs=%s\n" % ",".join(batch_sample_ids))

        suff_header = ["#CHROM", "POS", "REF", "Depth(CoveredSample)", "MappingQuality", "Readbases",
                       "ReadbasesQuality", "ReadPositionRank", "Strand"]
        OUT.write("%s\n" % "\t".join(suff_header))

        for i in range(region_number):

            batch_buffer = region_batch_buffers[i]
            position_number = len(batch_buffer.batch_heap)
            for j in range(position_number):
                OUT.write("%s\n" % batch_buffer.batch_heap[j])
    return


