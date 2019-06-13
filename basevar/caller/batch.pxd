"""Data structure for batch

Author: Shujia Huang
Date: 2019-06-05 10:28:22
"""
import cython

from basevar.io.fasta cimport FastaFile
from basevar.io.read cimport cAlignedRead


@cython.final
cdef class BatchElement:
    cdef public bytes ref_name
    cdef public long int ref_pos
    cdef public bytes ref_base
    cdef public bytes read_base  # Could be single base or an Indel
    cdef public int mapq  # mapping quality
    cdef public int base_qual  # A base quality for single base, set to be 0 if `read_base` is Indels
    cdef public int read_pos_rank
    cdef public char map_strand
    cdef public int base_type  # SIN(0), INS(1) or DEL(2)
    cdef void add_batch_element(self, BatchElement other)


@cython.final
cdef class BatchGenerator:
    cdef int CIGAR_M
    cdef int CIGAR_I
    cdef int CIGAR_D
    cdef int CIGAR_N
    cdef int CIGAR_S
    cdef int CIGAR_H
    cdef int CIGAR_P
    cdef int CIGAR_EQ
    cdef int CIGAR_X
    cdef int min_map_qual
    cdef int min_base_qual

    cdef dict batch_heap  # Main data will all be in ``batch_heap``

    cdef FastaFile ref_fa
    cdef bytes ref_name
    cdef bytes py_refseq
    cdef char* refseq
    cdef long int ref_seq_start
    cdef long int ref_seq_end
    cdef long int reg_start
    cdef long int reg_end
    cdef int qual_bin_size
    cdef object options

    cdef void create_batch_in_region(self, tuple region, cAlignedRead** read_start, cAlignedRead** read_end)
    cdef void add_batchelement_to_list(self, BatchElement be)
    cdef void _get_batch_from_single_read_in_region(self, cAlignedRead* read, long int start, long int end)
    cdef void _get_matchbase_from_read_segment(self, char* read_seq, char* read_qual, int mapq, char map_strand,
                                               int read_start, int read_offset, int ref_offset, int seglength,
                                               long int region_start, long int region_end)
