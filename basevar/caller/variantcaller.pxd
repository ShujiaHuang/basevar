"""Header for variantcaller.pyx
"""
cdef extern from "stdlib.h" nogil:
    int atoi(const char *str)
    void *calloc(size_t, size_t)
    void free(void *)

cdef extern from "string.h" nogil:
    char *strsep(char **string_ptr, const char *delimiter)

cdef class BatchInfo:
    cdef int size
    cdef bytes chrid
    cdef int position
    cdef bytes ref_base
    cdef int depth
    cdef int *mapqs
    cdef char *strands
    cdef char **sample_bases  # could be indel sequence
    cdef int *sample_base_quals
    cdef int *read_pos_rank
    cdef void destroy(self)

cdef bint variants_discovery(bytes chrid, list batchfiles, dict popgroup, float min_af,
                             int batch_count, cvg_file_handle, vcf_file_handle)