"""launch module

Author: Shujia Huang
Date: 2020-02-21 14:27:11
"""

from basevar.utils cimport BaseTypeCmdOptions
from basevar.datatype.strarray cimport StringArray


cdef class BaseTypeRunner:
    cdef StringArray align_files
    cdef StringArray sample_ids
    cdef dict pop_group

    cdef unsigned nCPU
    cdef basestring reference_file
    cdef basestring outvcf
    cdef basestring outcvg
    cdef list regions_for_each_process

    cdef BaseTypeCmdOptions options
    cdef void _set_cmd_options(self, args)

    cdef bint basevar_caller(self)

