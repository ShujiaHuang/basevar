"""Fast cython implementation of some windowing functions.
"""
from basevar.io.htslibWrapper cimport cAlignedRead


cdef class ReadArray:
    cdef cAlignedRead** array
    cdef cAlignedRead** window_start
    cdef cAlignedRead** window_end
    cdef int __size
    cdef int __capacity
    cdef int __longest_read

    cdef void append(self, cAlignedRead* value)
    cdef int get_size(self)
    cdef void set_window_pointers(self, int start, int end)
    cdef int get_length_of_longest_read(self)
    cdef int count_reads_covering_region(self, int start, int end)
    cdef void set_window_pointers_based_on_mate_pos(self, int start, int end)


cdef class BamReadBuffer:
    cdef char* sample

    cdef char* chrom
    cdef int chrom_id
    cdef long int start
    cdef long int end
    cdef int* filtered_read_counts_by_type
    cdef bint is_sorted
    cdef int window_start_base
    cdef int window_end_base
    cdef int max_reads
    cdef int min_map_qual
    cdef int min_base_qual
    cdef int trim_overlapping
    cdef int trim_soft_clipped
    cdef int verbosity

    cdef ReadArray reads
    cdef ReadArray bad_reads
    cdef cAlignedRead* last_read

    cdef void set_window_pointers(self, long long int start, long long int end, long long int refstart,
                                  long long int refend, char* refseq, int qual_bin_size)
    cdef void recompress_reads_in_current_window(self, long long int refstart, long long int refend,
                                                 char* refseq, int qual_bin_size, int compress_reads)
    cdef void add_read_to_buffer(self, cAlignedRead* the_read)
    cdef int count_improper_pairs(self)
    cdef int count_alignment_gaps(self)
    cdef int count_reads_covering_region(self, long long int start, long long int end)
    cdef void sort_reads(self)

    cdef ReadArray broken_mates
    cdef void sort_broken_mates(self)


