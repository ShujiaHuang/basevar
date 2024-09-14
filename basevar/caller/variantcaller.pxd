"""Header for variantcaller.pyx
"""
from basevar.datatype.strarray cimport StringArray

cdef bint variants_discovery(const char *chrid, const StringArray *batchfiles, dict popgroup, float min_af,
                             cvg_file_handle, vcf_file_handle)

