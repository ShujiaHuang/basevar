"""Data structure for batch

Author: Shujia Huang
Date: 2019-06-05 10:28:22
"""
from basevar.io.fasta cimport FastaFile
from basevar.io.read cimport cAlignedRead

cdef extern from "stdlib.h" nogil:
    void *malloc(size_t)
    void *calloc(size_t, size_t)
    void *realloc(void *, size_t)
    void free(void *)

cdef extern from "string.h":
    ctypedef int size_t
    char *strcpy(char *dest, char *src)
    char *strcat(char *dest, char *src)
    size_t strlen(char *s)


cdef class BatchInfo:
    # record the size of array: `*mapqs`==`*strands`==`**sample_bases` == `*sample_base_quals` == `*read_pos_rank`
    cdef int size

    cdef bytes chrid
    cdef long int position
    cdef bytes ref_base
    cdef int depth

    cdef char **sample_bases  # could be indel sequence
    cdef int *sample_base_quals
    cdef int *read_pos_rank
    cdef int *mapqs
    cdef char *strands
    cdef int *is_empty

    cdef void update_info_by_index(self, int index, bytes _target_chrom, long int _target_position, int mapq,
                                   char map_strand, char *read_base, int base_qual, int read_pos_rank)
    cdef basestring get_str(self)
    cdef void fill_empty(self)


ctypedef struct Element:
    pass


ctypedef struct BatchCigar:
    char *sample_bases_cigar
    char *sample_base_quals_cigar
    char *read_pos_rank_cigar
    char *mapqs_cigar
    char *strands_cigar


cdef class PositionBatchCigarArray:
    """A class just convert BatchInfo value into Element as a compress value to save memory!"""

    cdef bytes chrid
    cdef long int position
    cdef bytes ref_base
    cdef int depth

    # These element arrays are all in the same size
    cdef BatchCigar *array

    cdef int __size
    cdef int __capacity

    cdef void append(self, BatchInfo value)
    cdef int size(self)
    cdef BatchInfo batchcigar2batchinfo(self)


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

    cdef int sample_size
    cdef list batch_heap  # A list of BatchInfo and data will all be here
    cdef long int start_pos_in_batch_heap

    cdef FastaFile ref_fa
    cdef bytes ref_name
    cdef char *refseq
    cdef long int ref_seq_start
    cdef long int ref_seq_end
    cdef long int reg_start
    cdef long int reg_end
    cdef int *filtered_read_counts_by_type
    cdef object options

    cdef void create_batch_in_region(self, tuple region, cAlignedRead **read_start, cAlignedRead **read_end,
                                     int sample_index)
    cdef void get_batch_from_single_read_in_region(self, cAlignedRead *read, long int start, long int end,
                                                   int sample_index)
    cdef void _get_matchbase_from_read_segment(self, int sample_index, char*read_seq, char*read_qual, int mapq,
                                               char map_strand, int read_start, int read_offset, int ref_offset,
                                               int seglength, long int region_start, long int region_end)
