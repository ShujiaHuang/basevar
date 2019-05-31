from cpython cimport bool
cpdef list get_sample_names(list bamfiles, bool filename_has_samplename)
cdef list load_bamdata(list bamfiles, bytes chrom, int start, int end, char* refseq)
