"""Header for variantcaller.pyx
"""
cdef extern from "stdlib.h" nogil:
    int atoi(const char *str)
    void *calloc(size_t, size_t)
    void free(void *)

cdef extern from "string.h" nogil:
    char *strsep(char **string_ptr, const char *delimiter)

cdef bint variants_discovery(bytes chrid, list batchfiles, dict popgroup, float min_af,
                             int batch_count, cvg_file_handle, vcf_file_handle)