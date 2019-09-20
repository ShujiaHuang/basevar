"""Header for variantcaller.pyx
"""
cdef extern from "stdlib.h" nogil:
    int atoi(const char *str)
    void *calloc(size_t, size_t)
    void free(void *)

cdef extern from "string.h" nogil:
    char *strsep(char ** string_ptr, const char *delimiter)

from basevar.io.fasta cimport FastaFile

cdef bint variant_discovery_in_regions(FastaFile fa,
                                       list align_files,
                                       list regions,
                                       list samples,
                                       dict popgroup,
                                       basestring out_cvg_file_name,
                                       basestring out_vcf_file_name,
                                       object options)

