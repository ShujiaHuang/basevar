"""Fast cython implementation of some windowing functions.
"""
from basevar.utils cimport BaseTypeCmdOptions
from basevar.io.htslibWrapper cimport cAlignedRead
from basevar.datatype.genomeregion cimport GenomeRegion

cdef bint check_and_trim_read(cAlignedRead *the_read, cAlignedRead *the_last_read, int *filtered_read_counts_by_type,
                              int min_map_qual, bint trim_overlapping, bint trim_soft_clipped)

cdef class ReadArray:
    cdef cAlignedRead **array
    cdef cAlignedRead **window_start
    cdef cAlignedRead **window_end
    cdef unsigned long __size
    cdef unsigned long __capacity
    cdef unsigned __longest_read

    cdef void append(self, cAlignedRead *value)
    cdef unsigned long get_size(self)
    cdef void set_window_pointers(self, int start, int end)
    cdef unsigned get_length_of_longest_read(self)
    cdef unsigned long count_reads_covering_region(self, int start, int end)
    cdef void set_window_pointers_based_on_mate_pos(self, int start, int end)

cdef class BamReadBuffer:
    cdef char *sample

    cdef unsigned chrom_id
    cdef GenomeRegion region  # This GenomeRegion is 0-base.
    cdef int *filtered_read_counts_by_type

    cdef bint is_sorted
    cdef unsigned long window_start_base
    cdef unsigned long window_end_base

    cdef BaseTypeCmdOptions options

    cdef ReadArray reads
    cdef ReadArray bad_reads
    cdef cAlignedRead *last_read

    cdef void set_window_pointers(self, long int start, long int end, long int refstart,
                                  long int refend, char *refseq, int qual_bin_size)
    cdef void recompress_reads_in_current_window(self, long int refstart, long int refend, char *refseq,
                                                 int qual_bin_size, int compress_reads)
    cdef void add_read_to_buffer(self, cAlignedRead *the_read)
    cdef int count_improper_pairs(self)
    cdef int count_alignment_gaps(self)
    cdef int count_reads_covering_region(self, long int start, long int end)
    cdef void log_filter_summary(self)

    cdef ReadArray broken_mates

    # cdef void sort_reads(self)
    # cdef void sort_broken_mates(self)
