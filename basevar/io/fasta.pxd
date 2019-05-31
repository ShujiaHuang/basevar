"""IO for fasta file.

Author: Shujia Huang
Date: 2019-05-29
"""
cdef class FastaIndex:
    cdef long int n_targets
    cdef dict references
    cpdef dict target_name
    cpdef dict target_length


cdef class FastaFile:
    cdef bytes filename
    cdef object the_file
    cdef FastaIndex the_index

    cdef dict references
    cdef bytes cache
    cdef bytes cache_ref_name
    cdef long long int cache_start_pos
    cdef long long int cache_end_pos

    cpdef void close(self)
    cdef bytes get_character(self, bytes seq_name, long long int pos)
    cdef bytes get_sequence(self, bytes seq_name, long long int begin_pos, long long int end_pos)
    cdef void set_cache_sequence(self, bytes seq_name, long long int begin_pos, long long int end_pos)
