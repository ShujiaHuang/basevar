"""Package for batch file create

Author: Shujia Huang
Date: 2019-06-04 16:13:08
"""
from basevar.io.fasta cimport FastaFile
from basevar.utils cimport BaseTypeCmdOptions
from basevar.datatype.strarray cimport StringArray
from basevar.datatype.genomeregion cimport GenomeRegion

cdef void generate_batchfile(const GenomeRegion region,  # 1-base in GenomeRegion
                             const StringArray *batch_align_files,
                             const StringArray *batch_sample_ids,
                             FastaFile fa,
                             char *out_batch_file,
                             BaseTypeCmdOptions options)
