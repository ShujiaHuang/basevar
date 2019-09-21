"""Loading bam

Author: Shujia Huang
Date: 2019-06-03 00:28:50
"""
from basevar.io.htslibWrapper cimport Samfile
from basevar.caller.batch cimport BatchGenerator

cdef list get_sample_names(list bamfiles, bint filename_has_samplename)
cdef bint load_data_from_bamfile(Samfile bam_reader,
                                 bytes sample_id,
                                 bytes chrom,
                                 long int start, # 1-base
                                 long int end,   # 1-base
                                 BatchGenerator sample_batch_buffers,
                                 int sample_index,
                                 options)
