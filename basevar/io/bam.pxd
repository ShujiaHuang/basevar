"""Loading bam

Author: Shujia Huang
Date: 2019-06-03 00:28:50
"""
from basevar.io.htslibWrapper cimport Samfile
from basevar.caller.batch cimport BatchGenerator

from basevar.utils cimport BaseTypeCmdOptions
from basevar.datatype.strarray cimport StringArray
from basevar.datatype.genomeregion cimport GenomeRegion

cdef StringArray get_sample_names(StringArray bamfiles, bint filename_has_samplename)

cdef list load_bamdata(const StringArray *bamfiles,
                       const StringArray *samples,
                       const GenomeRegion region,
                       BaseTypeCmdOptions options)

cdef bint load_data_from_bamfile(Samfile bam_reader,
                                 char *sample_id,
                                 GenomeRegion region, # 1-base in GenomeRegion
                                 BatchGenerator sample_batch_buffers,
                                 unsigned long sample_index,
                                 options)

