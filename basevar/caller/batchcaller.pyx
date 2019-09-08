# cython: profile=True
"""
Package for parsing bamfile
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
import os
import time

from basevar.log import logger
from basevar.io.openfile import Open
from basevar.io.fasta cimport FastaFile
from basevar.caller.batch cimport BatchGenerator, BatchInfo
from basevar.io.read cimport check_and_trim_read
from basevar.io.htslibWrapper cimport Read_IsQCFail
from basevar.io.htslibWrapper cimport Samfile, ReadIterator, cAlignedRead

import sys
from basevar.io.bam cimport load_bamdata
from basevar.io.read cimport BamReadBuffer


cdef list create_batchfiles_in_regions(bytes chrom_name,
                                       list regions,
                                       long int region_boundary_start,
                                       long int region_boundary_end,
                                       list align_files,
                                       FastaFile fa,
                                       list samples,
                                       bytes outdir,
                                       object options):
    """
    ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
    ``samples``: The sample id of align_files
    ``fa``:
        # get sequence of chrom_name from reference fasta
        fa = self.ref_file_hd.fetch(chrid)

    """
    # store all the batch files
    cdef list batchfiles = []
    cdef int batchcount = options.batch_count
    cdef int part_num = len(align_files) / batchcount
    if part_num * batchcount < len(align_files):
        part_num += 1

    cdef bytes refseq_bytes = fa.get_sequence(chrom_name, region_boundary_start,
                                              region_boundary_end+1.0*options.r_len)
    cdef char* refseq = refseq_bytes

    cdef int m = 0
    cdef int i = 0
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
                           region_boundary_start,
                           region_boundary_end,
                           regions,
                           sub_align_files,
                           refseq,
                           fa,
                           options,
                           part_file_name,
                           batch_sample_ids)

        logger.info("Done for batchfile %s , %d seconds elapsed." % (
            part_file_name, time.time() - start_time))

    return batchfiles


cdef void generate_batchfile(bytes chrom_name,
                             long int bigstart,
                             long int bigend,
                             list regions,
                             list batch_align_files,
                             char* refseq,
                             FastaFile fa,
                             object options,
                             bytes out_batch_file,
                             list batch_sample_ids):

    """Loading bamfile and create a batchfile in ``regions``.

    Parameters:
        ``bigstart``: It's already 0-base position
        ``bigend``: It's already 0-base position
        ``regions``: The coordinate in regions is 1-base system.
    """
    # Just loading baminfo and not need to load header!
    cdef dict bamfiles = {s: f for s, f in zip(batch_sample_ids, batch_align_files)}
    cdef list read_buffers

    try:
        # load the whole mapping reads in [chrom_name, bigstart, bigend]
        read_buffers = load_bamdata(bamfiles, batch_sample_ids, chrom_name, bigstart, bigend, refseq, options)

    except Exception, e:
        logger.error("Exception in region %s:%s-%s. Error: %s" % (chrom_name, bigstart+1, bigend+1, e))
        sys.exit(1)

    if read_buffers is None or len(read_buffers) == 0:
        logger.info("Skipping region %s:%s-%s as it's empty." % (chrom_name, bigstart+1, bigend+1))
        return

    # take all the batch information for all samples in ``regions``
    cdef BamReadBuffer sample_read_buffer
    cdef BatchGenerator be_generator
    cdef list batch_buffers = []
    cdef int longest_read_size = 0
    cdef long int reg_start, reg_end
    for sample_read_buffer in read_buffers:

        if longest_read_size < sample_read_buffer.reads.get_length_of_longest_read():
            longest_read_size = sample_read_buffer.reads.get_length_of_longest_read()

        # init the batch generator by the BIG region, but the big region may not be use
        be_generator = BatchGenerator(chrom_name, bigstart, bigend, fa, options)

        # get batch information for each sample in regions
        for reg_start, reg_end in regions:
            reg_start -= 1  # set position to be 0-base
            reg_end -= 1  # set position to be 0-base
            be_generator.create_batch_in_region(
                (chrom_name, reg_start, reg_end),
                sample_read_buffer.reads.array,  # this start pointer will move automatically
                sample_read_buffer.reads.array + sample_read_buffer.reads.get_size()
            )

        batch_buffers.append(be_generator)

    # Todo: take care here, although these may not been called forever.
    if longest_read_size > options.r_len:
        options.r_len = longest_read_size

    output_batch_file(chrom_name, fa, batch_buffers, out_batch_file, batch_sample_ids, regions)
    return


# This function take much more memery than the ``generate_batchfile``, but I don't know why.
cdef void generate_batchfile_2(bytes chrom_name,
                               long int bigstart,
                               long int bigend,
                               list regions,
                               list batch_align_files,
                               char* refseq,
                               FastaFile fa,
                               object options,
                               bytes out_batch_file,
                               list batch_sample_ids):
    """Loading bamfile and create a batchfile in ``regions``.
    
    Parameters:
        ``bigstart``: It's already 0-base position
        ``bigend``: It's already 0-base position
        ``regions``: The coordinate in regions is 1-base system.
    """
    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead* the_read

    cdef BatchGenerator be_generator
    cdef list batch_buffers = []
    cdef int longest_read_size = 0
    cdef long int reg_start, reg_end

    cdef int mapq = options.mapq
    cdef bint trim_overlapping = options.trim_overlapping
    cdef bint trim_soft_clipped = options.trim_soft_clipped

    cdef bint read_ok
    cdef bint is_empty
    cdef bint is_overlap
    cdef int reg_size = len(regions)
    cdef int reg_index = 0
    cdef int i = 0

    _py_region = "%s:%s-%s" % (chrom_name, bigstart, bigend)
    cdef char* region = _py_region

    for sample, bamfile in zip(batch_sample_ids, batch_align_files):

        reader = Samfile(bamfile)
        reader.open("r", True)

        # init the batch generator by the BIG region, but the big region may not be use
        be_generator = BatchGenerator(chrom_name, bigstart, bigend, fa, options)

        reg_index = 0
        is_empty  = False

        try:
            reader_iter = reader.fetch(region)
        except Exception as e:
            logger.warning(e.message)
            logger.warning("No data could be retrieved for sample %s in file %s in "
                           "region %s" % (sample, reader.filename, region))
            is_empty = True

        while (not is_empty) and reader_iter.cnext():

            the_read = reader_iter.get(0, NULL)
            if the_read == NULL:
                continue

            read_ok = check_and_trim_read(the_read, NULL, be_generator.filtered_read_counts_by_type,
                                          mapq, trim_overlapping, trim_soft_clipped)
            if Read_IsQCFail(the_read):
                continue

            is_overlap = False
            for i in range(reg_index, reg_size):
                reg_start, reg_end = regions[i][0] -1, regions[i][1] - 1  # 0-base

                # still behind the region, do nothing but continue
                if the_read.end < reg_start:
                    continue

                # Break the loop when mapping start position is outside the region.
                if the_read.pos > reg_end:
                    break

                reg_index = i
                is_overlap = True
                break

            if is_overlap:
                reg_start, reg_end = regions[reg_index][0] -1, regions[reg_index][1] - 1  # 0-base
                be_generator.get_batch_from_single_read_in_region(the_read, reg_start, reg_end)
                if longest_read_size < the_read.end - the_read.pos:
                    longest_read_size = the_read.end - the_read.pos

        reader.close()
        batch_buffers.append(be_generator)

    # Todo: take care these codes, although they may not been called forever.
    if longest_read_size > options.r_len:
        options.r_len = longest_read_size

    output_batch_file(chrom_name, fa, batch_buffers, out_batch_file, batch_sample_ids, regions)
    return


cdef void output_batch_file(bytes chrom_name, FastaFile fa, list sample_batch_buffers, bytes out_batch_file,
                            list batch_sample_ids, list regions):

    cdef int sample_size = len(sample_batch_buffers)
    cdef BatchInfo batch_info = BatchInfo(chrom_name, sample_size)

    cdef BatchGenerator be_generator
    cdef bytes kk, tmp_bytes
    cdef int position, i
    with Open(out_batch_file, "wb", isbgz=True) if out_batch_file.endswith(".gz") else \
            open(out_batch_file, "w") as OUT:

        OUT.write("##fileformat=BaseVarBatchFile_v1.0\n")
        if batch_sample_ids:
            OUT.write("##SampleIDs=%s\n" % ",".join(batch_sample_ids))

        suff_header = ["#CHROM", "POS", "REF", "Depth(CoveredSample)", "MappingQuality", "Readbases",
                       "ReadbasesQuality", "ReadPositionRank", "Strand"]
        OUT.write("%s\n" % "\t".join(suff_header))

        for reg_start, reg_end in regions:
            for position in range(reg_start, reg_end+1):

                # `batch_info` will be update each loop here.
                # reset depth
                batch_info.depth = 0
                batch_info.position = position
                batch_info.ref_base = fa.get_character(chrom_name, position-1)  # 0-base system

                # set to be 0-base: position - 1
                kk = bytes("%s:%s" % (chrom_name, position-1))
                i = 0
                for i in range(sample_size):
                    be_generator = sample_batch_buffers[i]
                    if kk in be_generator.batch_heap:

                        batch_info.depth += 1
                        if be_generator.batch_heap[kk].base_type == 0:
                            # Single base
                            batch_info.sample_bases[i] = be_generator.batch_heap[kk].read_base
                        elif be_generator.batch_heap[kk].base_type == 1:
                            # insertion
                            tmp_bytes = "+"+be_generator.batch_heap[kk].read_base
                            batch_info.sample_bases[i] = tmp_bytes
                        elif be_generator.batch_heap[kk].base_type == 2:
                            # deletion
                            tmp_bytes = "-"+be_generator.batch_heap[kk].ref_base
                            batch_info.sample_bases[i] = tmp_bytes
                        else:
                            raise TypeError, ("Unknown base-type %s in 'output_batch_file'."
                                              "\n" % be_generator.batch_heap[kk].read_base)

                        batch_info.sample_base_quals[i] = be_generator.batch_heap[kk].base_qual
                        batch_info.strands[i] = be_generator.batch_heap[kk].map_strand
                        batch_info.mapqs[i] = be_generator.batch_heap[kk].mapq
                        batch_info.read_pos_rank[i] = be_generator.batch_heap[kk].read_pos_rank+1
                    else:
                        batch_info.sample_bases[i] = 'N'
                        batch_info.sample_base_quals[i] = 0
                        batch_info.strands[i] = '.'
                        batch_info.mapqs[i] = 0
                        batch_info.read_pos_rank[i] = 0

                OUT.write("%s\n" % batch_info)
    return
