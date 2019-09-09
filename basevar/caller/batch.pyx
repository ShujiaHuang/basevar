# cython: profile=True
"""
Author: Shujia Huang
Date: 2019-06-05 10:59:21
"""
import sys
import cython

from basevar.log import logger
from basevar.io.read cimport cAlignedRead
from basevar.io.htslibWrapper cimport Read_IsQCFail
from basevar.io.htslibWrapper cimport Read_IsReverse

cdef int SIN = 0   # Just a single base
cdef int INS = 1
cdef int DEL = 2

# the same with definition in read.pyx
cdef int LOW_QUAL_BASES = 0
cdef int UNMAPPED_READ = 1
cdef int MATE_UNMAPPED = 2
cdef int MATE_DISTANT = 3
cdef int SMALL_INSERT = 4
cdef int DUPLICATE = 5
cdef int LOW_MAP_QUAL = 6


cdef extern from "stdlib.h" nogil:
    void *calloc(size_t, size_t)
    void free(void *)


# A class to store batch information.
cdef class BatchInfo:

    def __cinit__(self, bytes chrid, int size):
        self.size = size
        self.chrid = chrid
        self.position = 0
        self.depth = 0
        self.strands = <char*>(calloc(size, sizeof(char)))
        self.sample_bases = <char**>(calloc(size, sizeof(char*)))
        self.sample_base_quals = <int*>(calloc(size, sizeof(int)))
        self.read_pos_rank = <int*>(calloc(size, sizeof(int)))
        self.mapqs = <int*>(calloc(size, sizeof(int)))

        # check initialization
        assert self.strands != NULL, "Could not allocate memory for self.strands in BatchInfo."
        assert self.sample_bases != NULL, "Could not allocate memory for self.sample_bases in BatchInfo."
        assert self.sample_base_quals != NULL, "Could not allocate memory for self.sample_base_quals in BatchInfo."
        assert self.read_pos_rank != NULL, "Could not allocate memory for self.read_pos_rank in BaseInfo."
        assert self.mapqs != NULL, "Could not allocate memory for self.mapqs in BaseInfo."

    # def __str__(self):
    #     return self.get_str()
    #
    # cdef basestring get_str(self):
    #     """
    #     __str__ is called when you do str(BatchInfo) or print BatchInfo, and will return a short string, which
    #     describing the BatchInfo.
    #     """
    #     cdef list sample_bases      = []
    #     cdef list sample_base_quals = []
    #     cdef list read_pos_rank     = []
    #     cdef list mapqs             = []
    #     cdef list strands           = []
    #
    #     cdef int i = 0
    #     if self.depth > 0:
    #         for i in range(self.size):
    #             mapqs.append(str(self.mapqs[i]))
    #             sample_bases.append(str.upper(self.sample_bases[i]))  # charater to upper()
    #             sample_base_quals.append(str(self.sample_base_quals[i]))
    #             read_pos_rank.append(str(self.read_pos_rank[i]))
    #             strands.append(chr(self.strands[i]))
    #
    #         return "\t".join([
    #             self.chrid,
    #             str(self.position),
    #             self.ref_base,
    #             str(self.depth),
    #             ",".join(mapqs),
    #             ",".join(sample_bases),
    #             ",".join(sample_base_quals),
    #             ",".join(read_pos_rank),
    #             ",".join(strands)
    #         ])
    #     else:
    #         return "\t".join(map(str, [self.chrid, self.position, self.ref_base, self.depth, ".\t.\t.\t.\t."]))

    def __dealloc__(self):
        self.destroy()

    cdef void destroy(self):
        """Free memory"""

        self.depth = 0
        if self.strands != NULL:
            free(self.strands)

        cdef int i = 0
        if self.sample_bases != NULL:
            free(self.sample_bases)

        if self.sample_base_quals != NULL:
            free(self.sample_base_quals)

        if self.read_pos_rank != NULL:
            free(self.read_pos_rank)

        if self.mapqs != NULL:
            free(self.mapqs)

        return


cdef class BatchElement(object):
    """
    Class to encapsulate information for all position. The basic idea is to tread all position as batch.
    """
    def __cinit__(self, bytes ref_name, long int ref_pos, bytes ref_base, bytes read_base,
                  int mapq, int base_qual, int read_pos_rank, char map_strand):

        self.ref_name = ref_name
        self.ref_pos = ref_pos
        self.ref_base = ref_base
        self.read_base = read_base
        self.base_qual = base_qual
        self.read_pos_rank = read_pos_rank
        self.mapq = mapq
        self.map_strand = map_strand

        if len(self.ref_base) == len(read_base):
            self.base_type = SIN
        elif len(self.read_base) == 0:
            self.base_type = DEL
        elif len(self.ref_base) == 0:
            self.base_type = INS
        else:
            raise TypeError, "Unknown base-type %s\n" % read_base

    def __str__(self):
        """
        __str__ is called when you do str(BatchElement), and will return a short string, which
        describing the BatchElement
        """
        str_info = "%s:%s" % (self.ref_name, self.ref_pos)
        return str_info

    # Todo: Update this function later!
    cdef void add_batch_element(self, BatchElement other):
        pass


cdef class BatchGenerator(object):
    """
    A class to generate batch informattion from a bunch of reads.
    """
    def __cinit__(self, bytes ref_name, long int reg_start, long int reg_end, FastaFile ref_fa,
                  options):
        """
        Constructor. Create a storage place for batchfile, and store the values of some flags which
        are used in the pysam CIGAR information.
        """
        self.CIGAR_M  = 0 # Match
        self.CIGAR_I  = 1 # Insertion
        self.CIGAR_D  = 2 # Deletion
        self.CIGAR_N  = 3 # Skipped region from reference
        self.CIGAR_S  = 4 # Soft clipping. Sequence is present in read
        self.CIGAR_H  = 5 # Hard clipping. Sequence is not present in read
        self.CIGAR_P  = 6 # Padding. Used for padded alignment
        self.CIGAR_EQ = 7 # Alignment match; sequence match
        self.CIGAR_X  = 8 # Alignment match; sequence mismatch

        self.batch_heap = {} # List of batch position
        self.options    = options

        self.ref_fa        = ref_fa
        self.ref_name      = ref_name
        self.reg_start     = reg_start  # start position of region, and region[1] must be 0-base coordinate system
        self.reg_end       = reg_end  # end position of region, and region[2] must be 0-base coordinate system

        self.ref_seq_start = max(0, self.reg_start-200)
        self.ref_seq_end   = min(self.reg_end+200, self.ref_fa.references[self.ref_name].seq_length-1)
        self.py_refseq     = self.ref_fa.get_sequence(self.ref_name, self.ref_seq_start, self.ref_seq_end) # Cache this
        self.refseq        = self.py_refseq

        # the same definition with class ``BamReadBuffer`` in read.pyx
        self.filtered_read_counts_by_type = <int*>(calloc(7, sizeof(int)))
        if options.filter_duplicates == 0:
            self.filtered_read_counts_by_type[DUPLICATE] = -1

        if options.filter_reads_with_unmapped_mates == 0:
            self.filtered_read_counts_by_type[MATE_UNMAPPED] = -1

        if options.filter_reads_with_distant_mates == 0:
            self.filtered_read_counts_by_type[MATE_DISTANT] = -1

        if options.filterReadPairsWithSmallInserts == 0:
            self.filtered_read_counts_by_type[SMALL_INSERT] = -1

    def __dealloc__(self):
        """Clean up memory.
        """
        if self.filtered_read_counts_by_type != NULL:
            free(self.filtered_read_counts_by_type)

    cdef void add_batchelement_to_list(self, BatchElement be):
        """Check if the element is already in the heap. if not, add it. """

        cdef BatchElement the_be = self.batch_heap.get(be, None)
        if the_be is None:
            self.batch_heap[str(be)] = be
        else:
            # do nothing, which means we just get the first base from the first alignment read
            pass
            # the_be.add_batch_element(be)

    cdef void create_batch_in_region(self, tuple region, cAlignedRead** read_start, cAlignedRead** read_end):
        """Fetch batch information in a specific region."""
        cdef bytes chrom = region[0]
        cdef long int start = region[1]  # must be 0-base coordinate system
        cdef long int end = region[2]  # must be 0-base coordinate system

        assert chrom == self.ref_name, "Error match chromosome (%s != %s)" % (chrom, self.ref_name)

        if start < self.ref_seq_start:
            logger.error("Start position (%s) is outside the reference region (%s)" % (
                start, "%s:%s-%s" % (self.ref_name, self.ref_seq_start, self.ref_seq_end)))
            sys.exit(1)

        if end > self.ref_seq_end:
            logger.error("End position (%s) is outside the reference region (%s)" % (
                end, "%s:%s-%s" % (self.ref_name, self.ref_seq_start, self.ref_seq_end)))
            sys.exit(1)

        cdef int read_num = 0
        while read_start != read_end:

            if Read_IsQCFail(read_start[0]):
                read_start += 1 # QC fail read move to the next one
                continue

            # still behind the region, do nothing but continue
            if read_start[0].end < start:
                read_start += 1
                continue

            # Break the loop when mapping start position is outside the region.
            if read_start[0].pos > end:
                break

            # get batch information here!
            self.get_batch_from_single_read_in_region(read_start[0], start, end)

            read_num += 1 # how many reads in this regions
            read_start += 1  # move to the next read

        return

    cdef void get_batch_from_single_read_in_region(self, cAlignedRead* read, long int start, long int end):
        """Check a single read for batch. 
        Batch positions are flagged by the CIGAR string. Pysam reports the CIGAR string information as a list 
        of tuples, where each tuple is a pair, and the first element gives the type of feature (match, insertion 
        or deletion), and the second element gives the number of nucleotides associated. 
        For example, [(0, 1), (1, 2), (0, 1)] is a 1 base match, a 2 base insertion, and a 1 base match.
        """
        cdef long int read_start_pos = read.pos

        cdef int ref_offset    = 0
        cdef int read_offset   = 0
        cdef int mapq          = read.mapq
        cdef int cigar_length  = read.cigar_len
        cdef char* read_seq    = read.seq
        cdef char* read_qual   = read.qual
        cdef bytes insert_seq  = None
        cdef bytes deleted_seq = None
        cdef char map_strand   = '-' if Read_IsReverse(read) else '+'

        cdef int cigar_flag  = 0
        cdef int cigar_index = 0
        cdef int length = 0
        cdef long int ref_pos = 0

        for cigar_index in range(cigar_length):
            cigar_flag = read.cigar_ops[2*cigar_index]
            length = read.cigar_ops[(2*cigar_index)+1]

            # An insertion take us further along the read, but not the reference
            if cigar_flag == self.CIGAR_I:

                if cigar_index > 0 and read.cigar_ops[(2*cigar_index)-2] == self.CIGAR_M:
                    pass
                elif cigar_index < cigar_length - 1 and read.cigar_ops[(2*cigar_index)+2] == self.CIGAR_M:
                    pass
                else:
                    read_offset += length
                    continue

                insert_seq = read_seq[read_offset:read_offset+length]
                ref_pos = read_start_pos + ref_offset - 1
                # do not use insertion with Ns in them
                if (start <= ref_pos <= end) and (insert_seq.count("N") == 0):

                    self.add_batchelement_to_list(BatchElement(
                        self.ref_name, ref_pos, "", insert_seq, mapq, 0, read_offset, map_strand
                    ))

                read_offset += length

            # A deletion take us further along the reference, but not the read
            elif cigar_flag == self.CIGAR_D:
                if cigar_index > 0 and read.cigar_ops[(2*cigar_index)-2] == self.CIGAR_M:
                    pass
                elif cigar_index < cigar_length -1 and read.cigar_ops[(2*cigar_index)+2] == self.CIGAR_M:
                    pass
                else:
                    ref_offset += length
                    continue

                deleted_seq = self.ref_fa.get_sequence(self.ref_name, read_start_pos+ref_offset,
                                                       read_start_pos+ref_offset+length)
                ref_pos = read_start_pos + ref_offset - 1
                # do not use deletion with Ns in them
                if (start <= ref_pos <= end) and (deleted_seq.upper().count("N") == 0):
                    self.add_batchelement_to_list(BatchElement(
                        self.ref_name, ref_pos, deleted_seq, "", mapq, 0, read_offset, map_strand
                    ))

                ref_offset += length

            # A match take us further along the reference and the read
            elif cigar_flag == self.CIGAR_M or cigar_flag == self.CIGAR_EQ or cigar_flag == self.CIGAR_X:

                self._get_matchbase_from_read_segment(read_seq, read_qual, mapq, map_strand, read_start_pos,
                                                      read_offset, ref_offset, length, start, end)

                ref_offset += length
                read_offset += length

            elif cigar_flag == self.CIGAR_N:
                ref_offset += length

            elif cigar_flag == self.CIGAR_S:
                read_offset += length
                # We have to move back read position when there is a soft-clipping at the
                # beginning of reads
                if cigar_index == 0:
                    ref_offset += length

            # Hard clipping. Sequence is not present in read.
            elif cigar_flag == self.CIGAR_H:
                continue

            # Padding. We do nothing here.
            elif cigar_flag == self.CIGAR_P:
                continue

            # Other kinds of flag. Just don't care about them.
            else:
                continue

        return

    cdef void _get_matchbase_from_read_segment(self, char* read_seq, char* read_qual, int mapq, char map_strand,
                                               int read_start, int read_offset, int ref_offset, int seglength,
                                               long int region_start, long int region_end):
        """Get BatchElement from a particular CIGAR_M segment in read.
        
        Parameter
        =========
            ``read_seq``: The complete sequence of bases in the read
            ``read_qual``: The complete sequence of quality scores in the read
            ``read_start``: The starting position, in the reference sequence, of the read
            ``read_offset``: int
                If we're not starting from the beginning of the read, then how far along 
                the read sequence to start.
            ``ref_offset``: int
                If the read contains indels, we need to offset our reference position accordingly, 
                by this much.
        """
        # ignore all the positions which are not in [region_start, region_end]
        if (read_start + ref_offset - 1 > region_end) or (read_start + ref_offset + seglength < region_start):
            return

        cdef int base_qual = 0
        cdef int read_index = 0
        cdef int ref_index = 0
        cdef bytes read_char
        cdef bytes ref_char

        cdef int index = 0
        cdef long int ref_pos = 0
        for index in range(seglength):

            ref_pos = read_start + ref_offset + index
            if ref_pos < region_start:
                continue

            if ref_pos > region_end:
                break

            read_index = read_offset + index
            read_char = read_seq[read_index]
            base_qual = read_qual[read_index]

            ref_index = ref_pos - self.ref_seq_start
            ref_char = self.refseq[ref_index]

            self.add_batchelement_to_list(BatchElement(
                self.ref_name, ref_pos, ref_char, read_char, mapq, base_qual, read_index, map_strand
            ))

        return
