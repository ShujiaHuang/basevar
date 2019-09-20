"""This is a Process module for BaseType by BAM/CRAM

Author: Shujia Huang
Date: 2019-06-10 09:19:19
"""
from basevar.io.fasta cimport FastaFile

cdef class BaseVarProcess:
    cdef list samples
    cdef list align_files

    cdef FastaFile fa_file_hd
    cdef list regions
    cdef dict dict_regions
    cdef dict popgroup

    cdef basestring out_vcf_file
    cdef basestring out_cvg_file
    # cdef bytes cache_dir

    cdef object options

    # cdef void run_variant_discovery_by_batch(self)
    cdef void run_variant_discovery_in_regions(self)
