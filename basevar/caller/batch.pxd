"""Data structure for batch

Author: Shujia Huang
Date: 2019-06-05 10:28:22
"""
from basevar.io.fasta cimport FastaFile
from basevar.io.read cimport cAlignedRead

cdef extern from "stdlib.h" nogil:
    void *calloc(size_t, size_t)
    void free(void *)

cdef extern from "string.h":
    ctypedef int size_t
    char *strcpy(char *dest, char *src)
    char *strcat(char *dest, char *src)
    size_t strlen(char *s)

ctypedef struct Element:
    int n    # number
    char *b   # charater or string

cdef class BatchInfo:
    # record the size of array: `*mapqs`==`*strands`==`**sample_bases` == `*sample_base_quals` == `*read_pos_rank`
    cdef int size

    cdef bytes chrid
    cdef long int position
    cdef bytes ref_base
    cdef int depth
    cdef int *mapqs
    cdef char *strands
    cdef char **sample_bases  # could be indel sequence
    cdef int *sample_base_quals
    cdef int *read_pos_rank
    cdef int *is_empty

    cdef void update_info_by_index(self, int index, bytes _target_chrom, long int _target_position, int mapq,
                                   char map_strand, char *read_base, int base_qual, int read_pos_rank)
    cdef basestring get_str(self)
    cdef void fill_empty(self)

cdef class StaticBatchElement:
    cdef Element sample_bases
    cdef Element sample_base_quals
    cdef Element read_pos_ranks
    cdef Element strands
    cdef Element mapqs

    cdef int __size
    cdef int __capacity

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
