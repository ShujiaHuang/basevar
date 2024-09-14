"""Data structure for batch

Author: Shujia Huang
Date: 2019-06-05 10:28:22
"""
from basevar.io.fasta cimport FastaFile
from basevar.io.read cimport cAlignedRead
from basevar.utils cimport BaseTypeCmdOptions
from basevar.datatype.genomeregion cimport GenomeRegion
from basevar.datatype.dynamicstring cimport dstring


cdef class BatchInfo:
    # record the size of array: `*mapqs`==`*strands`==`**sample_bases` == `*sample_base_quals` == `*read_pos_rank`
    cdef size_t size
    cdef size_t __capacity

    cdef char *chrid
    cdef long int position
    cdef char *ref_base
    cdef long int depth

    cdef char **sample_bases  # could be indel sequence
    cdef int *sample_base_quals
    cdef int *read_pos_rank
    cdef int *mapqs
    cdef char *strands
    cdef int *is_empty

    cdef void set_empty(self)
    cdef void set_size(self, size_t size)
    cdef void clear(self)
    cdef void update_info_by_index(self, int index, char *_target_chrom, unsigned long _target_position, int mapq,
                                   char map_strand, char *read_base, int base_qual, int read_pos_rank)
    cdef dstring get_str(self)

# ``_Cigar_string``, ``_Cigar_char`` and ``_Cigar_int`` are all for matching the type
# in ``BatchInfo``
ctypedef struct _CigarString:
    int n
    char *b

ctypedef struct _CigarInt:
    int n
    int b

ctypedef struct _CigarChar:
    int n
    char b

ctypedef struct _CigarStringArray:
    int size
    _CigarString *data

ctypedef struct _CigarCharArray:
    int size
    _CigarChar *data

ctypedef struct _CigarIntArray:
    int size
    _CigarInt *data

ctypedef struct BatchCigar:
    # string(or char**) type
    _CigarStringArray sample_bases_cigar

    # int type
    _CigarIntArray sample_base_quals_cigar
    _CigarIntArray read_pos_rank_cigar
    _CigarIntArray mapqs_cigar

    # single char type
    _CigarCharArray strands_cigar

cdef class PositionBatchCigarArray:
    """A class convert BatchInfo value into Element as a compress value to save memory!"""
    cdef char *chrid
    cdef unsigned long position
    cdef char *ref_base

    cdef BatchCigar *array

    cdef unsigned long __size
    cdef unsigned long __capacity
    cdef unsigned long __sample_number
    cdef unsigned long __depth

    cdef unsigned long size(self)
    cdef void append(self, BatchInfo value)
    cdef BatchInfo convert_position_batch_cigar_array_to_batchinfo(self) # convert the whole array into one BatchInfo

    cdef _CigarStringArray _compress_string(self, char **data, int size)
    cdef _CigarCharArray _compress_char(self, char *data, int size)
    cdef _CigarIntArray _compress_int(self, int *data, int size)
    cdef BatchCigar _BatchInfo2BatchCigar(self, BatchInfo value)


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
    cdef long int start_pos_in_batch_heap
    cdef list batch_heap

    cdef GenomeRegion region # 1-base system
    cdef FastaFile ref_fa
    cdef long int ref_seq_start
    cdef long int ref_seq_end
    cdef int *filtered_read_counts_by_type
    cdef BaseTypeCmdOptions options

    cdef void create_batch_in_region(self, const GenomeRegion region, cAlignedRead **read_start, cAlignedRead **read_end,
                                     int sample_index)
    cdef void get_batch_from_single_read_in_region(self, cAlignedRead *read, long int start, long int end,
                                                   int sample_index)
    cdef void _get_matchbase_from_read_segment(self, int sample_index,
                                               const char *read_seq,
                                               const char *read_qual,
                                               int mapq,
                                               char map_strand,
                                               int read_start,
                                               int read_offset,
                                               int ref_offset,
                                               int seglength,
                                               long int region_start,
                                               long int region_end)
