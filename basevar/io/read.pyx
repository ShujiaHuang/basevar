# cython: profile=True
"""Fast cython implementation of some windowing functions.
"""
from basevar.log import logger
from basevar.io.htslibWrapper cimport cAlignedRead
from basevar.io.htslibWrapper cimport destroy_read
from basevar.io.htslibWrapper cimport compress_read
from basevar.io.htslibWrapper cimport uncompress_read
from basevar.io.htslibWrapper cimport Read_IsReverse
from basevar.io.htslibWrapper cimport Read_IsPaired
from basevar.io.htslibWrapper cimport Read_IsProperPair
from basevar.io.htslibWrapper cimport Read_IsDuplicate
from basevar.io.htslibWrapper cimport Read_IsUnmapped
from basevar.io.htslibWrapper cimport Read_MateIsUnmapped
from basevar.io.htslibWrapper cimport Read_MateIsReverse
from basevar.io.htslibWrapper cimport Read_IsSecondaryAlignment
from basevar.io.htslibWrapper cimport Read_SetQCFail
from basevar.io.htslibWrapper cimport Read_IsCompressed


cdef int LOW_QUAL_BASES = 0
cdef int UNMAPPED_READ = 1
cdef int MATE_UNMAPPED = 2
cdef int MATE_DISTANT = 3
cdef int SMALL_INSERT = 4
cdef int DUPLICATE = 5
cdef int LOW_MAP_QUAL = 6


cdef extern from "stdlib.h":
    void *malloc(size_t)
    void *calloc(size_t, size_t)
    void *realloc(void *, size_t)
    void free(void *)
    void qsort(void*, size_t, size_t, int(*)(void*, void*))


cdef extern from "math.h":
    int abs(int)


# @cython.profile(False)
# cdef inline int read_pos_comp(const void* x, const void* y) nogil:
#     """
#     Comparison function for use in qsort, to sort reads by their start positions and
#     then end positions.
#     """
#     cdef cAlignedRead** read_one = <cAlignedRead**> (x)
#     cdef cAlignedRead** read_two = <cAlignedRead**> (y)
#
#     return read_one[0].pos - read_two[0].pos
#
#
# @cython.profile(False)
# cdef int read_mate_pos_comp(const void* x, const void* y):
#     """
#     Comparison function for use in qsort, to sort reads by their mate positions, as long as
#     the mates are on the same chromosome.
#     """
#     cdef cAlignedRead** read_one = <cAlignedRead**> (x)
#     cdef cAlignedRead** read_two = <cAlignedRead**> (y)
#
#     # Sorting is broken for reads with mates on different chromosomes
#     assert read_one[0].mate_chrom_id == read_two[0].mate_chrom_id
#     return read_one[0].mate_pos - read_two[0].mate_pos


cdef class ReadArray:
    """Simple structure to wrap a raw C array, with some bounds checking.
    """
    def __cinit__(self, int size):
        """Allocate an array of size 'size', with initial values 'init'.
        """
        self.array = <cAlignedRead**> (malloc(size * sizeof(cAlignedRead*)))
        assert self.array != NULL, "Could not allocate memory for ReadArray"

        self.__size = 0  # We don't put anything in here yet, just allocate memory
        self.__capacity = size
        self.__longest_read = 0
        self.window_start = NULL
        self.window_end = NULL

        # Always initialise to NULL
        cdef int index = 0
        for index in range(size):
            self.array[index] = NULL

    def __dealloc__(self):
        """
        Free memory
        """
        cdef int index = 0
        if self.array != NULL:
            for index in range(self.__size):
                if self.array[index] != NULL:
                    destroy_read(self.array[index])
                    self.array[index] = NULL

            free(self.array)

    cdef int get_size(self):
        """
        Return the size of the array
        """
        return self.__size

    cdef void append(self, cAlignedRead* value):
        """Append a new value to the array, re-allocating if necessary.
        """
        cdef cAlignedRead** temp = NULL

        if self.__size == self.__capacity:
            temp = <cAlignedRead**>(realloc(self.array, 2 * sizeof(cAlignedRead*) * self.__capacity))

            if temp == NULL:
                raise StandardError, "Could not re-allocate ReadArray"
            else:
                self.array = temp
                self.__capacity *= 2

        self.array[self.__size] = value
        self.__size += 1

        cdef int read_length = value.end - value.pos
        if read_length > self.__longest_read:
            self.__longest_read = read_length

    cdef int count_reads_covering_region(self, int start, int end):
        """
        Return the number of reads which overlap this region, where 'end' is not
        included, and 'start' is.
        """
        cdef cAlignedRead **r_start
        cdef cAlignedRead **r_end
        cdef int first_overlap_start = -1
        cdef int start_pos_of_reads = -1
        cdef int end_pos_of_reads = -1

        # Set pointers for good reads
        if self.__size == 0:
            return 0
        else:
            first_overlap_start = max(1, start - self.__longest_read)
            start_pos_of_reads = bisect_reads_left(self.array, first_overlap_start, self.__size)
            end_pos_of_reads = bisect_reads_left(self.array, end, self.__size)

            while start_pos_of_reads < self.__size and self.array[start_pos_of_reads].end <= start:
                start_pos_of_reads += 1

            r_start = self.array + start_pos_of_reads
            r_end = min(self.array + end_pos_of_reads, self.array + self.__size)

            if start_pos_of_reads > end_pos_of_reads:
                logger.error("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" % (
                    start, end, start_pos_of_reads, end_pos_of_reads))
                logger.error("There are %s reads here." % (self.__size))
                raise StandardError, "This should never happen. Read start pointer > read end pointer!!"

            return r_end - r_start

    cdef void set_window_pointers(self, int start, int end):
        """
        Set the pointers 'windowStart' and 'windowEnd' to
        point to the relevant first and last +1 reads as specified by
        the co-ordinates.
        """
        cdef int first_overlap_start = -1
        cdef int start_pos_of_reads = -1
        cdef int end_pos_of_reads = -1

        # Set pointers for good reads
        if self.__size == 0:
            self.window_start = self.array
            self.window_end = self.array
        else:
            first_overlap_start = max(1, start - self.__longest_read)
            start_pos_of_reads = bisect_reads_left(self.array, first_overlap_start, self.__size)
            end_pos_of_reads = bisect_reads_left(self.array, end, self.__size)

            while start_pos_of_reads < self.__size and self.array[start_pos_of_reads].end <= start:
                start_pos_of_reads += 1

            self.window_start = self.array + start_pos_of_reads
            self.window_end = min(self.array + end_pos_of_reads, self.array + self.__size)

            if start_pos_of_reads > end_pos_of_reads:
                logger.info("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" % (
                start, end, start_pos_of_reads, end_pos_of_reads))
                logger.info("There are %s reads here." % (self.__size))
                raise StandardError, "This should never happen. Read start pointer > read end pointer!!"

    cdef void set_window_pointers_based_on_mate_pos(self, int start, int end):
        """
        Set the pointers 'window_start' and 'window_end' to
        point to the relevant first and last +1 reads as specified by
        the co-ordinates of the mates of the reads in this array.
        """
        cdef int first_overlap_start = -1
        cdef int start_pos_of_reads = -1
        cdef int end_pos_of_reads = -1

        # Set pointers for good reads
        if self.__size == 0:
            self.window_start = self.array
            self.window_end = self.array
        else:
            first_overlap_start = max(1, start - self.__longest_read)
            start_pos_of_reads = bisect_reads_left(self.array, first_overlap_start, self.__size, 1)
            end_pos_of_reads = bisect_reads_left(self.array, end, self.__size, 1)

            self.window_start = self.array + start_pos_of_reads
            self.window_end = min(self.array + end_pos_of_reads, self.array + self.__size)

            if start_pos_of_reads > end_pos_of_reads:
                logger.info("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" % (
                start, end, start_pos_of_reads, end_pos_of_reads))
                logger.info("There are %s reads here." % (self.__size))
                raise StandardError, "This should never happen. Read start pointer > read end pointer!!"

    cdef int get_length_of_longest_read(self):
        """Return the longest read size
        """
        return self.__longest_read


cdef int bisect_reads_left(cAlignedRead** reads, int test_pos, int n_reads, int test_mate_pos=0):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = n_reads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if not test_mate_pos:
            if reads[mid].pos < test_pos:
                low = mid + 1
            else:
                high = mid
        else:
            if reads[mid].mate_pos < test_pos:
                low = mid + 1
            else:
                high = mid

    return low


cdef bint check_and_trim_read(cAlignedRead* the_read, cAlignedRead* the_last_read, int* filtered_read_counts_by_type,
                             int min_map_qual, bint trim_overlapping, bint trim_soft_clipped):
    """
    Performs various quality checks on the read, and trims read (i.e. set q-scores to zero). Returns
    true if read is ok, and false otherwise.
    """
    if Read_IsSecondaryAlignment(the_read):
        Read_SetQCFail(the_read)
        return False

    if the_read.mapq < min_map_qual:
        filtered_read_counts_by_type[LOW_MAP_QUAL] += 1
        Read_SetQCFail(the_read)
        return False

    # Remove unmapped reads
    if Read_IsUnmapped(the_read):
        filtered_read_counts_by_type[UNMAPPED_READ] += 1
        Read_SetQCFail(the_read)
        return False

    # Remove broken pairs, i.e. pairs where the mate is mapped to a different chromosome or the
    # mate is unmapped
    if filtered_read_counts_by_type[MATE_UNMAPPED] != -1:

        if Read_IsPaired(the_read) and Read_MateIsUnmapped(the_read):
            filtered_read_counts_by_type[MATE_UNMAPPED] += 1
            return False

    if filtered_read_counts_by_type[MATE_DISTANT] != -1:

        if Read_IsPaired(the_read) and (
                the_read.chrom_id != the_read.mate_chrom_id or (not Read_IsProperPair(the_read))):

            filtered_read_counts_by_type[MATE_DISTANT] += 1
            return False

    # If the insert size is < read length, then we almost certainly have adapter contamination, in which
    # case we should skip these reads, as they may be mapped to the wrong location in the genome.
    if filtered_read_counts_by_type[SMALL_INSERT] != -1:

        if Read_IsPaired(the_read) and (the_read.insert_size != 0 and abs(the_read.insert_size) < the_read.r_len):
            filtered_read_counts_by_type[SMALL_INSERT] += 1
            Read_SetQCFail(the_read)
            return False

    # Check if this read is actually a duplicate.
    # Todo: store library tag and check.
    if filtered_read_counts_by_type[DUPLICATE] != -1:

        if Read_IsDuplicate(the_read):
            filtered_read_counts_by_type[DUPLICATE] += 1
            Read_SetQCFail(the_read)
            return False

        elif the_last_read != NULL:
            if the_read.pos == the_last_read.pos and the_read.r_len == the_last_read.r_len:

                # For paired reads, check mate's position
                if Read_IsPaired(the_read):

                    if the_last_read.mate_pos == the_read.mate_pos:
                        filtered_read_counts_by_type[DUPLICATE] += 1
                        Read_SetQCFail(the_read)
                        return False

                # For single reads, just check pos and length of reads
                else:
                    filtered_read_counts_by_type[DUPLICATE] += 1
                    Read_SetQCFail(the_read)
                    return False

    cdef int abs_ins = abs(the_read.insert_size)

    # Trim overlapping part of forward read, in pairs where the read length is greater than the insert size
    # N.B Insert size is from start of forward read to end of reverse read, i.e. fragment size. This is done to
    # remove duplicate information, which gives systematic errors when pcr errors have occured in library prep.
    cdef int i = 0
    if trim_overlapping and (
            Read_IsPaired(the_read) and abs_ins > 0 and (not Read_IsReverse(the_read)) and
            Read_MateIsReverse(the_read) and abs_ins < 2 * the_read.r_len):

        for i in range(1, min(the_read.r_len, (2 * the_read.r_len - the_read.insert_size) + 1)):
            the_read.qual[the_read.r_len - i] = 0

    cdef int cigar_index = 0
    cdef int cigar_op = -1
    cdef int cigar_len = -1

    cdef int index = 0
    cdef int j = 0
    if trim_soft_clipped:
        # Check for soft-clipping (present in BWA reads, for example, but not Stampy). Soft-clipped
        # sequences should be set to QUAL = 0, as they may include contamination by adapters etc.

        for cigar_index in range(the_read.cigar_len):

            cigar_op = the_read.cigar_ops[2 * cigar_index]
            cigar_len = the_read.cigar_ops[2 * cigar_index + 1]

            # Skip good sequence. 0 is match. 1 is insertion.
            if cigar_op == 0 or cigar_op == 1:
                index += the_read.cigar_ops[2 * cigar_index + 1]

            # 4 is soft-clipping
            elif cigar_op == 4:
                # Set quals to zero across all sequence flagged as soft-clipped
                for j in range(the_read.cigar_ops[2 * cigar_index + 1]):
                    the_read.qual[index] = 0
                    index += 1

    return True

cdef class BamReadBuffer:
    """
    Utility class for bufffering reads from a single BAM file, so we only make a single pass
    through the data in each BAM in the loop through windows.
    """
    def __cinit__(self, char* chrom, long int start, long int end, options):
        """
        Constructor.

        ``start``: 0-base
        ``end``: 0-base
        """
        cdef int initial_size = max(100, ((end - start)/options.r_len))
        self.is_sorted = True

        self.reads = ReadArray(initial_size)
        self.bad_reads = ReadArray(initial_size)
        self.broken_mates = ReadArray(initial_size)
        self.filtered_read_counts_by_type = <int*>(calloc(7, sizeof(int)))
        self.chrom = chrom
        self.start = start
        self.end = end

        # self.max_reads = options.max_reads
        self.min_map_qual = options.mapq
        self.trim_overlapping = options.trim_overlapping
        self.trim_soft_clipped = options.trim_soft_clipped
        self.verbosity = options.verbosity

        self.last_read = NULL

        if options.filter_duplicates == 0:
            self.filtered_read_counts_by_type[DUPLICATE] = -1

        if options.filter_reads_with_unmapped_mates == 0:
            self.filtered_read_counts_by_type[MATE_UNMAPPED] = -1

        if options.filter_reads_with_distant_mates == 0:
            self.filtered_read_counts_by_type[MATE_DISTANT] = -1

        if options.filter_read_pairs_with_small_inserts == 0:
            self.filtered_read_counts_by_type[SMALL_INSERT] = -1

    def __dealloc__(self):
        """Clean up memory.
        """
        if self.filtered_read_counts_by_type != NULL:
            free(self.filtered_read_counts_by_type)

    cdef void log_filter_summary(self):
        """Useful debug information about which reads have been filtered out.
        """
        if self.verbosity >= 3:
            region = "%s:%s-%s" % (self.chrom, self.start+1, self.end+1)
            logger.debug("Sample %s has %s good reads in %s" % (self.sample, self.reads.get_size(), region))
            logger.debug("Sample %s has %s bad reads in %s" % (self.sample, self.bad_reads.get_size(), region))
            logger.debug("Sample %s has %s broken mates %s" % (self.sample, self.broken_mates.get_size(), region))
            logger.debug("|-- low map quality reads = %s" % (self.filtered_read_counts_by_type[LOW_MAP_QUAL]))
            logger.debug("|-- low qual reads = %s" % (self.filtered_read_counts_by_type[LOW_QUAL_BASES]))
            logger.debug("|-- un-mapped reads = %s" % (self.filtered_read_counts_by_type[UNMAPPED_READ]))
            logger.debug("|-- reads with unmapped mates = %s" % (self.filtered_read_counts_by_type[MATE_UNMAPPED]))
            logger.debug("|-- reads with distant mates = %s" % (self.filtered_read_counts_by_type[MATE_DISTANT]))
            logger.debug("|-- reads pairs with small inserts = %s" % (self.filtered_read_counts_by_type[SMALL_INSERT]))
            logger.debug("|__ duplicate reads = %s\n" % (self.filtered_read_counts_by_type[DUPLICATE]))

            if self.trim_overlapping == 1:
                logger.debug("Overlapping segments of read pairs were clipped")
            else:
                logger.debug("Overlapping segments of read pairs were not clipped")

            logger.debug("Overhanging bits of reads were not clipped")

    cdef void add_read_to_buffer(self, cAlignedRead *the_read):
        """Add a new read to the buffer, making sure to re-allocate memory when necessary.
        """
        cdef int read_ok = 0
        cdef int min_good_bases_this_read = 0
        cdef int read_start = -1
        cdef int read_end = -1
        cdef int read_length = -1

        if the_read == NULL:
            return
        else:

            # TODO: Check that this works for duplicates when first read goes into bad reads pile...
            if self.last_read != NULL:
                read_ok = check_and_trim_read(the_read, self.last_read, self.filtered_read_counts_by_type,
                                              self.min_map_qual, self.trim_overlapping,
                                              self.trim_soft_clipped)

                if self.last_read.pos > the_read.pos:
                    self.is_sorted = False
            else:
                read_ok = check_and_trim_read(the_read, NULL, self.filtered_read_counts_by_type, self.min_map_qual,
                                              self.trim_overlapping, self.trim_soft_clipped)

            # ignore read which the same mapping position
            if self.reads.get_size() > 0 and self.last_read.pos == the_read.pos:
                return

            # self.last_read = the_read
            # Put read into bad array
            if not read_ok:
                self.bad_reads.append(the_read)

            # Put read into good array
            else:
                self.last_read = the_read
                self.reads.append(the_read)

    cdef int count_alignment_gaps(self):
        """
        Count and return the number of indels seen
        by the mapper in all good and bad reads.
        """
        cdef cAlignedRead** start = self.reads.window_start
        cdef cAlignedRead** end = self.reads.window_end
        cdef cAlignedRead** b_start = self.bad_reads.window_start
        cdef cAlignedRead** b_end = self.bad_reads.window_end

        cdef int n_gaps = 0
        cdef int i = 0

        while start != end:
            for i in range(start[0].cigar_len):
                if 1 <= start[0].cigar_ops[2 * i] <= 4:
                    n_gaps += 1

            start += 1

        while b_start != b_end:
            for i in range(b_start[0].cigar_len):
                if 1 <= b_start[0].cigar_ops[2 * i] <= 4:
                    n_gaps += 1

            b_start += 1

        return n_gaps

    cdef int count_improper_pairs(self):
        """
        Count and return the number of reads (Good and bad) that
        are members of improper pairs.
        """
        cdef cAlignedRead** start = self.reads.window_start
        cdef cAlignedRead** end = self.reads.window_end
        cdef cAlignedRead** b_start = self.bad_reads.window_start
        cdef cAlignedRead** b_end = self.bad_reads.window_end

        cdef int n_improper = 0

        while start != end:
            if not Read_IsProperPair(start[0]):
                n_improper += 1
            start += 1

        while b_start != b_end:
            if not Read_IsProperPair(b_start[0]):
                n_improper += 1
            b_start += 1

        return n_improper

    cdef int count_reads_covering_region(self, long long int start, long long int end):
        """
        Return the number of 'good' reads covering this region.
        """
        return self.reads.count_reads_covering_region(start, end)

    cdef void set_window_pointers(self, long long int start, long long int end, long long int refstart,
                                  long long int refend, char* refseq, int qual_bin_size):
        """
        Set the window_start and window_end pointers to point to the first
        and last+1 reads covering this window.
        """
        self.reads.set_window_pointers(start, end)
        self.bad_reads.set_window_pointers(start, end)
        self.broken_mates.set_window_pointers_based_on_mate_pos(start, end)

        cdef cAlignedRead** the_start = NULL
        cdef cAlignedRead** the_end = NULL

        the_start = self.reads.window_start
        the_end = self.reads.window_end

        # If the reads in this window are compressed, then uncompress them
        if the_start != the_end and Read_IsCompressed(the_start[0]):

            while the_start != the_end:
                uncompress_read(the_start[0], refseq, refstart, refend, qual_bin_size)
                the_start += 1

            the_start = self.bad_reads.window_start
            the_end = self.bad_reads.window_end

            while the_start != the_end:
                uncompress_read(the_start[0], refseq, refstart, refend, qual_bin_size)
                the_start += 1

            the_start = self.broken_mates.window_start
            the_end = self.broken_mates.window_end

            while the_start != the_end:
                uncompress_read(the_start[0], refseq, refstart, refend, qual_bin_size)
                the_start += 1

    cdef void recompress_reads_in_current_window(self, long long int refstart, long long int refend,
                                                 char* refseq, int qual_bin_size, int is_compress_reads):
        """
        Set the windowStart and windowEnd pointers to point to the first
        and last+1 reads covering this window.
        """
        cdef cAlignedRead** the_start = NULL
        cdef cAlignedRead** the_end = NULL

        # If the reads in this window are compressed, then uncompress them
        the_start = self.reads.window_start
        the_end = self.reads.window_end
        while the_start != the_end:

            # Compresss reads is required
            if is_compress_reads:
                compress_read(the_start[0], refseq, refstart, refend, qual_bin_size)

            # Otherwise, just clear the hash
            if the_start[0].hash != NULL:
                free(the_start[0].hash)
                the_start[0].hash = NULL

            the_start += 1

        # For bad_reads
        the_start = self.bad_reads.window_start
        the_end = self.bad_reads.window_end
        while the_start != the_end:

            # Compresss reads is required
            if is_compress_reads:
                compress_read(the_start[0], refseq, refstart, refend, qual_bin_size)

            # Otherwise, just clear the hash
            if the_start[0].hash != NULL:
                free(the_start[0].hash)
                the_start[0].hash = NULL

            the_start += 1

        the_start = self.broken_mates.window_start
        the_end = self.broken_mates.window_end
        while the_start != the_end:

            # Compresss reads is required
            if is_compress_reads:
                compress_read(the_start[0], refseq, refstart, refend, qual_bin_size)

            # Otherwise, just clear the hash
            if the_start[0].hash != NULL:
                free(the_start[0].hash)
                the_start[0].hash = NULL

            the_start += 1

    # cdef void sort_reads(self):
    #     """
    #     Sort the contents of the reads array by position
    #     """
    #     if not self.is_sorted:
    #         if self.reads.get_size() > 0:
    #             qsort(self.reads.array, self.reads.get_size(), sizeof(cAlignedRead**), read_pos_comp)
    #
    #         if self.bad_reads.get_size() > 0:
    #             qsort(self.bad_reads.array, self.bad_reads.get_size(), sizeof(cAlignedRead**),
    #                   read_pos_comp)
    #
    #         self.is_sorted = True
    #
    # cdef void sort_broken_mates(self):
    #     """
    #     Sort the contents of the brokenMates array by the co-ordinates of their
    #     mates.
    #     """
    #     qsort(self.broken_mates.array, self.broken_mates.get_size(), sizeof(cAlignedRead**),
    #           read_mate_pos_comp)
