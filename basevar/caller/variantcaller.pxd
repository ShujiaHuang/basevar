"""Header for variantcaller.pyx
"""
cdef extern from "stdlib.h" nogil:
    int atoi(const char *str)
    void *calloc(size_t, size_t)
    void free(void *)

cdef extern from "string.h" nogil:
    char *strsep(char ** string_ptr, const char *delimiter)


from basevar.io.fasta cimport FastaFile

cdef bint variant_discovery_in_regions(FastaFile fa, list align_files, list regions, list samples, dict popgroup,
                                       bytes outdir, object options, out_cvg_file, out_vcf_file)

# cdef bint variants_discovery(bytes chrid, list batchfiles, dict popgroup, float min_af,
#                              int batch_count, cvg_file_handle, vcf_file_handle)
