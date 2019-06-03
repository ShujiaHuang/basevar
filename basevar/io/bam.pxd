from cpython cimport bool
cpdef list get_sample_names(list bamfiles, bool filename_has_samplename)
cdef list load_bamdata(dict bam_objs, list samples, bytes chrom, long long int start, long long int end,
                       char* refseq, options)
