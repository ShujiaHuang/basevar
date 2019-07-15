"""
Package for parsing bamfile
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
import sys
import os

from basevar.log import logger
from basevar.io.openfile import Open
from basevar.io.fasta cimport FastaFile
from basevar.io.bam cimport load_bamdata
from basevar.caller.batch cimport BatchGenerator
from basevar.io.read cimport BamReadBuffer
from basevar.io.htslibWrapper cimport compress_read
from basevar.io.htslibWrapper cimport Samfile, ReadIterator, cAlignedRead


cdef list create_batchfiles_in_regions(bytes chrom_name,
                                       list regions,
                                       long int region_boundary_start,
                                       long int region_boundary_end,
                                       list align_files,
                                       FastaFile fa,
                                       list samples,
                                       bytes outdir,
                                       object options,
                                       bint is_smart_rerun):
    """
    ``samples``: The sample id of align_files
    ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
    ``fa``:
        # get sequence of chrom_name from reference fasta
        fa = self.ref_file_hd.fetch(chrid)

    justbase: bool
                Just outout base in batch file
    """
    # store all the batch files
    cdef list batchfiles = []
    cdef int batchcount = options.batch_count
    cdef int part_num = len(align_files) / batchcount
    if part_num * batchcount < len(align_files):
        part_num += 1

    cdef bytes refseq_bytes = fa.get_sequence(
        chrom_name, region_boundary_start, region_boundary_end + 5 * options.r_len)
    cdef char* refseq = refseq_bytes

    cdef int m = 0
    cdef int i = 0
    for i in range(0, len(align_files), batchcount):
        # Create a batch of temp files which we call them batchfiles for variant discovery

        m += 1
        part_file_name = "basevar.%s.%d_%d.batch.gz" % (
            "%s.%s.%s" % (chrom_name, region_boundary_start+1, region_boundary_end+1), m, part_num)
        part_file_name = os.path.join(outdir, part_file_name)  # Join Path could fix different OS

        # store the name of batchfiles into a list.
        batchfiles.append(part_file_name)
        if is_smart_rerun and os.path.isfile(part_file_name):
            # ``part_file_name`` is exists We don't have to create it again if setting `smartrerun`
            logger.info("%s already exists, we don't have to create it again, "
                        "when you set `smartrerun`" % part_file_name)
            continue
        else:
            logger.info("Creating batchfile %s\n" % part_file_name)

        # One batch of alignment files, the size and order are the same between ``sub_align_files`` and
        # ``batch_sample_ids``
        sub_align_files = align_files[i:i + batchcount]
        batch_sample_ids = None
        if samples:
            batch_sample_ids = samples[i:i + batchcount]

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

    return batchfiles


cdef void generate_batchfile_old(bytes chrom_name,
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
        be_generator = BatchGenerator((chrom_name, bigstart, bigend), fa, options.mapq,
                                      options.min_base_qual, options)

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

    # Todo: take care these codes, although they may not been called forever.
    if longest_read_size > options.r_len:
        options.r_len = longest_read_size

    output_batch_file(chrom_name, fa, batch_buffers, out_batch_file,
                      batch_sample_ids, regions)
    return


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
    cdef Samfile reader
    cdef ReadIterator reader_iter
    cdef cAlignedRead* the_read

    cdef bint is_compress_read = options.is_compress_read
    cdef int qual_bin_size = options.qual_bin_size
    cdef int max_read_thd = options.max_reads

    cdef int total_reads = 0
    cdef BamReadBuffer sample_read_buffer

    cdef BatchGenerator be_generator
    cdef list batch_buffers = []
    cdef int longest_read_size = 0
    cdef long int reg_start, reg_end

    cdef bint is_empty = False
    region = "%s:%s-%s" % (chrom_name, bigstart, bigend)
    for sample, bamfile in zip(batch_sample_ids, batch_align_files):

        reader = Samfile(bamfile)
        reader.open("r", True)

        # Need to lock the thread here when sharing BAM files
        if reader.lock is not None:
            reader.lock.acquire()

        # set initial size for BamReadBuffer
        sample_read_buffer = BamReadBuffer(chrom_name, bigstart, bigend, options)
        sample_read_buffer.sample = sample

        is_empty = False
        try:
            reader_iter = reader.fetch(region)
        except Exception as e:
            logger.warning(e.message)
            logger.warning("No data could be retrieved for sample %s in file %s in "
                           "region %s" % (sample, reader.filename, region))

            is_empty = True

        while (not is_empty) and reader_iter.cnext():

            the_read = reader_iter.get(0, NULL)

            if is_compress_read:
                compress_read(the_read, refseq, bigstart, bigend, qual_bin_size)

            sample_read_buffer.add_read_to_buffer(the_read)

            total_reads += 1
            if total_reads > max_read_thd:
                logger.error("Too many reads (%s) in region %s. Quitting now. Either reduce --buffer-size or "
                             "increase --max_reads." % (total_reads, region))

                reader.close()
                sys.exit(1)

        reader.close()  # close bamfile

        # Need to release thread lock here when sharing BAM files
        if reader.lock is not None:
            reader.lock.release()

        # population_read_buffers.append(sample_read_buffer)
        if longest_read_size < sample_read_buffer.reads.get_length_of_longest_read():
            longest_read_size = sample_read_buffer.reads.get_length_of_longest_read()

        # init the batch generator by the BIG region, but the big region may not be use
        be_generator = BatchGenerator((chrom_name, bigstart, bigend), fa, options.mapq,
                                      options.min_base_qual, options)

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

    # Todo: take care these codes, although they may not been called forever.
    if longest_read_size > options.r_len:
        options.r_len = longest_read_size

    output_batch_file(chrom_name, fa, batch_buffers, out_batch_file,
                      batch_sample_ids, regions)
    return


cdef void output_batch_file(bytes chrom_name, FastaFile fa, list batch_buffers, bytes out_batch_file,
                            list batch_sample_ids, list regions):
    cdef int depth = 0
    cdef BatchGenerator be_generator
    cdef list sample_bases, sample_base_quals, strands, mapqs, read_pos_rank
    cdef long int position
    cdef bytes ref_base
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

                # clear data
                depth = 0
                sample_bases = []
                sample_base_quals = []
                strands = []
                mapqs = []
                read_pos_rank = []

                # set to be 0-base: position - 1
                kk = "%s:%s" % (chrom_name, position-1)
                ref_base = fa.get_character(chrom_name, position-1)
                for be_generator in batch_buffers:

                    if kk in be_generator.batch_heap:
                        depth += 1
                        if be_generator.batch_heap[kk].base_type == 0:
                            # Single base
                            sample_bases.append(be_generator.batch_heap[kk].read_base)
                        elif be_generator.batch_heap[kk].base_type == 1:
                            # insertion
                            sample_bases.append("+"+be_generator.batch_heap[kk].read_base)
                        elif be_generator.batch_heap[kk].base_type == 2:
                            # deletion
                            sample_bases.append("-"+be_generator.batch_heap[kk].ref_base)
                        else:
                            raise TypeError, ("Unknown base-type %s in 'output_batch_file'."
                                              "\n" % be_generator.batch_heap[kk].read_base)

                        sample_base_quals.append(be_generator.batch_heap[kk].base_qual)
                        strands.append(chr(be_generator.batch_heap[kk].map_strand))
                        mapqs.append(be_generator.batch_heap[kk].mapq)
                        read_pos_rank.append(be_generator.batch_heap[kk].read_pos_rank+1)
                    else:
                        sample_bases.append("N")
                        sample_base_quals.append(0)
                        strands.append(".")
                        mapqs.append(0)
                        read_pos_rank.append(0)

                # OUT.write("%s\n" % "\t".join([
                #     chrom_name,
                #     str(position),
                #     ref_base,
                #     str(depth),
                #     ",".join(map(str, mapqs)),
                #     ",".join(sample_bases),
                #     ",".join(map(str, sample_base_quals)),
                #     ",".join(map(str, read_pos_rank)),
                #     ",".join(strands)
                # ]))

                OUT.write("%s\n" % "\t".join([
                    chrom_name,
                    str(position),
                    ref_base,
                    str(depth),
                    ",".join(map(str, mapqs)) if depth else ".",
                    ",".join(sample_bases) if depth else ".",
                    ",".join(map(str, sample_base_quals)) if depth else ".",
                    ",".join(map(str, read_pos_rank)) if depth else ".",
                    ",".join(strands) if depth else "."
                ]))

    return
