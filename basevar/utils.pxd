"""header for utils.pyx"""
from basevar.datatype.strarray cimport StringArray

ctypedef struct BaseTypeCmdOptions:
    # Not the whole list of basetype parameters.
    # Just some args options which will pass into c-type functions
    unsigned mapq                 # minimum mapping quality
    unsigned min_base_qual        # minimum base quality
    unsigned batch_count          # The number of samples in one batch file
    unsigned r_len                # The longest read length
    unsigned smartrerun

    unsigned long max_reads       # Maximum reads cover in a window
    unsigned is_compress_read     # If set to 1, then reads will be compressed and reduce memory usage
    unsigned qual_bin_size
    unsigned trim_overlapping
    unsigned trim_soft_clipped
    unsigned filter_duplicates
    unsigned filter_reads_with_unmapped_mates
    unsigned filter_reads_with_distant_mates
    unsigned filter_read_pairs_with_small_inserts

    unsigned verbosity
    double min_af                 # minimum AF threshold

cdef dict load_popgroup_info(const StringArray *samples, in_popgroup_file)
cdef long int c_max(long int x, long int y)
cdef long int c_min(long int x, long int y)
cdef void fast_merge_files(list temp_file_names, basestring final_file_name, bint is_del_raw_file)
cdef list generate_regions_by_process_num(list regions, int process_num, bint convert_to_2d)
