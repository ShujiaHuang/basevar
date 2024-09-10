"""Header for vcf_concat.pyx
"""
cdef extern from "stdlib.h" nogil:
    void *calloc(size_t, size_t)
    void free(void *)

cdef extern from "bcftools/bcftools.h":
    pass

cdef extern from "bcftools/vcfconcat.c":
    int main_vcfconcat(int argc, char *argv[])

cdef int call_vcfconcat(list argv)
