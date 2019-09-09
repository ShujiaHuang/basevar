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
    cdef dict popgroup

    cdef bytes out_vcf_file
    cdef bytes out_cvg_file
    cdef bytes cache_dir

    cdef object options
