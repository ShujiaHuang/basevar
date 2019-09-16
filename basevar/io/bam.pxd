"""Loading bam

Author: Shujia Huang
Date: 2019-06-03 00:28:50
"""

from basevar.caller.batch cimport BatchGenerator

cdef list get_sample_names(list bamfiles, bint filename_has_samplename)
cdef list load_bamdata(dict bamfiles, list samples, bytes chrom, long int start, long int end,
                       char* refseq, options)
cdef bint load_data_from_bamfile(dict bamfiles,
                                 list samples,
                                 bytes chrom,
                                 long int start, # 1-base
                                 long int end,   # 1-base
                                 BatchGenerator sample_batch_buffers,
                                 options)
