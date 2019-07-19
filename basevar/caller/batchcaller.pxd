"""Package for batch file create

Author: Shujia Huang
Date: 2019-06-04 16:13:08
"""
from basevar.io.fasta cimport FastaFile


cdef list create_batchfiles_in_regions(bytes chrom_name,
                                       list regions,
                                       long int region_boundary_start,
                                       long int region_boundary_end,
                                       list align_files,
                                       FastaFile fa,
                                       list sample_ids,
                                       bytes outdir,
                                       object options)