"""This is a Process module for BaseType by BAM/CRAM

Author: Shujia Huang
Date: 2019-06-10 09:19:19
"""
from basevar.utils cimport BaseTypeCmdOptions
from basevar.io.fasta cimport FastaFile
from basevar.datatype.strarray cimport StringArray
from basevar.datatype.genomeregion cimport GenomeRegionArray

cdef class BaseVarProcess:
    cdef StringArray align_files
    cdef StringArray samples

    cdef FastaFile fa_file_hd
    cdef GenomeRegionArray regions
    cdef dict pop_group

    cdef char *cache_dir
    cdef basestring out_vcf_file
    cdef basestring out_cvg_file

    cdef BaseTypeCmdOptions options
    cdef void run_variant_discovery_by_batchfiles(self)
