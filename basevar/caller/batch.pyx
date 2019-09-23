# cython: profile=True
"""
Author: Shujia Huang
Date: 2019-06-05 10:59:21
"""
import sys

from basevar.log import logger
from basevar.io.read cimport cAlignedRead
from basevar.io.htslibWrapper cimport Read_IsQCFail
from basevar.io.htslibWrapper cimport Read_IsReverse

cdef int SIN = 0  # Just a single base
cdef int INS = 1
cdef int DEL = 2

# the same with definition in read.pyx
cdef int LOW_QUAL_BASES = 0
cdef int UNMAPPED_READ = 1
cdef int MATE_UNMAPPED = 2
cdef int MATE_DISTANT = 3
cdef int SMALL_INSERT = 4
cdef int DUPLICATE = 5
cdef int LOW_MAP_QUAL = 6

# A class to store batch information.
cdef class BatchInfo:
    def __cinit__(self, bytes chrid, long int position=0, bytes ref_base=b'N', int size=0):
        self.size = size  # default size is as the same as capacity
        self.__capacity = size
        self.chrid = chrid
        self.position = position
        self.ref_base = ref_base

        self.depth = 0
        self.strands = <char*> (calloc(self.__capacity, sizeof(char)))
        self.sample_bases = <char**> (calloc(self.__capacity, sizeof(char*)))
        self.sample_base_quals = <int*> (calloc(self.__capacity, sizeof(int)))
        self.mapqs = <int*> (calloc(self.__capacity, sizeof(int)))
        self.read_pos_rank = <int*> (calloc(self.__capacity, sizeof(int)))
        self.is_empty = <int*> (calloc(self.__capacity, sizeof(int)))

        # check initialization
        assert self.sample_bases != NULL, "Could not allocate memory for self.sample_bases in BatchInfo."
        assert self.sample_base_quals != NULL, "Could not allocate memory for self.sample_base_quals in BatchInfo."
        assert self.read_pos_rank != NULL, "Could not allocate memory for self.read_pos_rank in BaseInfo."
        assert self.mapqs != NULL, "Could not allocate memory for self.mapqs in BaseInfo."
        assert self.strands != NULL, "Could not allocate memory for self.strands in BatchInfo."
        assert self.is_empty != NULL, "Could not allocate memory for self.is_empty in BaseInfo."

        # initialization and set empty mark for all the element, 1=>empty, 0=> not empty
        cdef int i
        for i in range(self.__capacity):
            self.is_empty[i] = 1
            self.mapqs[i] = 0
            self.strands[i] = '.'
            self.sample_bases[i] = 'N'
            self.sample_base_quals[i] = 0
            self.read_pos_rank[i] = 0

    def __str__(self):
        """
        __str__ is called when you do str(BatchInfo) or print BatchInfo, and will return a short string, which
        describing the BatchInfo.
        """
        return self.get_str()

    cdef void set_size(self, int size):
        """call this function to adjust the self.size, if the real size is not equal to capacity"""
        if size > self.__capacity:
            logger.error("The size (%d) is larger than the capacity(%d) in ``BatchInfo``." % (
                size, self.__capacity))
            sys.exit(1)

        self.size = size
        return

    cdef void set_empty(self):
        self.depth = 0

        cdef int i
        for i in range(self.__capacity):
            self.is_empty[i] = 1
            self.mapqs[i] = 0
            self.strands[i] = '.'
            self.sample_bases[i] = 'N'
            self.sample_base_quals[i] = 0
            self.read_pos_rank[i] = 0

        return

    cdef basestring get_str(self):

        cdef list sample_bases = []
        cdef list sample_base_quals = []
        cdef list read_pos_rank = []
        cdef list mapqs = []
        cdef list strands = []

        cdef int i = 0
        if self.depth > 0:
            for i in range(self.size):
                sample_bases.append(self.sample_bases[i])  # Note: sample_bases must all be upper
                sample_base_quals.append(str(self.sample_base_quals[i]))
                read_pos_rank.append(str(self.read_pos_rank[i]))
                strands.append(chr(self.strands[i]))
                mapqs.append(str(self.mapqs[i]))

            return "\t".join([self.chrid,
                              str(self.position),
                              self.ref_base,
                              str(self.depth),
                              ",".join(mapqs),
                              ",".join(sample_bases),
                              ",".join(sample_base_quals),
                              ",".join(read_pos_rank),
                              ",".join(strands)])
        else:
            return "\t".join(map(str, [self.chrid, self.position, self.ref_base, self.depth, ".\t.\t.\t.\t."]))

    cdef void update_info_by_index(self, int index, bytes _target_chrom, long int _target_position, int mapq,
                                   char map_strand, char *read_base, int base_qual, int read_pos_rank):
        """Update information"""
        assert _target_chrom == self.chrid and _target_position == self.position, \
            "Error! Chromosome(%s, %s) or position(%s, %s) not exactly match " % (
                self.chrid, _target_chrom, self.position, _target_position)

        if read_base[0] != '-' and read_base[0] != '+':
            self.depth += 1

        if self.is_empty[index]:
            self.is_empty[index] = 0  # not empty

        self.mapqs[index] = mapq
        self.strands[index] = map_strand
        self.sample_base_quals[index] = base_qual
        self.read_pos_rank[index] = read_pos_rank

        cdef int base_size = strlen(read_base)
        self.sample_bases[index] = <char*> (calloc(base_size, sizeof(char)))
        strcpy(self.sample_bases[index], read_base)

        return

    cdef void clear(self):
        # free memory
        self.__dealloc__()
        return

    def __dealloc__(self):
        """Free memory"""
        self.depth = 0
        if self.strands != NULL:
            free(self.strands)

        if self.sample_bases != NULL:
            free(self.sample_bases)

        if self.sample_base_quals != NULL:
            free(self.sample_base_quals)

        if self.read_pos_rank != NULL:
            free(self.read_pos_rank)

        if self.mapqs != NULL:
            free(self.mapqs)

        return

# compress the ``BatchInfo`` of ``BatchGenerator``
cdef class PositionBatchCigarArray:
    """A class for Element record of each position."""
    def __cinit__(self, bytes chrid, long int position, bytes ref_base, int array_size):
        """Allocate an array of size 'size', with initial values 'init'.
        """
        # base information at specific position!
        self.chrid = chrid
        self.position = position
        self.ref_base = ref_base

        self.array = <BatchCigar*>(malloc(array_size * sizeof(BatchCigar)))
        if self.array == NULL:
            raise StandardError, "Could not allocate memory for PositionBatchCigarArray"

        self.__size = 0  # We don't put anything in here yet
        self.__capacity = array_size
        self.__sample_number = 0  # The number of samples which have been store in this array
        self.__depth = 0

    def __dealloc__(self):
        """
        Free memory
        """
        cdef BatchCigar batch_cigar
        cdef int i = 0
        if self.array != NULL:
            for i in range(self.__size):
                batch_cigar = self.array[i]
                free(batch_cigar.sample_bases_cigar.data)
                free(batch_cigar.sample_base_quals_cigar.data)
                free(batch_cigar.read_pos_rank_cigar.data)
                free(batch_cigar.mapqs_cigar.data)
                free(batch_cigar.strands_cigar.data)

            free(self.array)

    cdef int size(self):
        return self.__size

    cdef void append(self, BatchInfo value):

        if value.chrid != self.chrid or value.position != self.position:
            raise StandardError, (
                    "Error in ``PositionBatchCigarArray`` when call append()! "
                    "Chromosome(%s, %s) or position(%s, %s) not exactly match " % (
                        self.chrid, value.chrid, self.position, value.position)
            )

        cdef BatchCigar *temp = NULL
        if self.__size == self.__capacity:
            temp = <BatchCigar*> (realloc(self.array, 2 * sizeof(BatchCigar) * self.__capacity))

            if temp == NULL:
                raise StandardError, "Could not re-allocate PositionBatchCigarArray!!"
            else:
                self.array = temp
                self.__capacity *= 2

        self.__depth += value.depth  # store coverage
        self.__sample_number += value.size

        self.array[self.__size] = self._BatchInfo2BatchCigar(value)
        self.__size += 1  # increase size

    cdef BatchCigar _BatchInfo2BatchCigar(self, BatchInfo value):
        # Set BatchCigar to compress BatchInfo and save memory.
        cdef BatchCigar batch_cigar

        batch_cigar.mapqs_cigar = self._compress_int(value.mapqs, value.size)
        batch_cigar.sample_bases_cigar = self._compress_string(value.sample_bases, value.size)

        batch_cigar.sample_base_quals_cigar = self._compress_int(value.sample_base_quals, value.size)
        batch_cigar.read_pos_rank_cigar = self._compress_int(value.read_pos_rank, value.size)
        batch_cigar.strands_cigar = self._compress_char(value.strands, value.size)

        return batch_cigar

    cdef BatchInfo convert_position_batch_cigar_array_to_batchinfo(self):

        cdef BatchInfo batch_info = BatchInfo(self.chrid, self.position, self.ref_base, self.__sample_number)
        batch_info.depth = self.__depth

        cdef BatchCigar batch_cigar
        cdef _CigarString cs
        cdef int i = 0, j = 0, _ = 0

        # index in batch_info
        cdef int m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0
        cdef int base_size = 0  # just for sample_bases

        cdef int total_base_array_size = 0
        cdef int total_qual_array_size = 0
        cdef int total_mapqs_array_size = 0
        cdef int total_pos_rank_array_size = 0
        cdef int total_strand_array_size = 0

        for i in range(self.__size):
            batch_cigar = self.array[i]

            # These total_* are just for debug
            total_base_array_size += batch_cigar.sample_bases_cigar.size
            total_qual_array_size += batch_cigar.sample_base_quals_cigar.size
            total_mapqs_array_size += batch_cigar.mapqs_cigar.size
            total_pos_rank_array_size += batch_cigar.read_pos_rank_cigar.size
            total_strand_array_size += batch_cigar.strands_cigar.size

            # set ``sample_bases``
            for j in range(batch_cigar.sample_bases_cigar.size):
                # print ">> ", self.chrid, self.position, self.ref_base, batch_cigar.sample_bases_cigar.data[j].n, batch_cigar.sample_bases_cigar.data[j].b
                if strcmp(batch_cigar.sample_bases_cigar.data[j].b, "N") != 0:

                    for _ in range(batch_cigar.sample_bases_cigar.data[j].n):
                        # copy char* type, must use strcpy().
                        base_size = strlen(batch_cigar.sample_bases_cigar.data[j].b)
                        batch_info.sample_bases[m1] = <char*>(calloc(base_size, sizeof(char)))
                        strcpy(batch_info.sample_bases[m1], batch_cigar.sample_bases_cigar.data[j].b)

                        m1 += 1
                else:
                    m1 += batch_cigar.sample_bases_cigar.data[j].n

            # set ``sample_base_quals``
            for j in range(batch_cigar.sample_base_quals_cigar.size):

                if batch_cigar.sample_base_quals_cigar.data[j].b != 0:
                    for _ in range(batch_cigar.sample_base_quals_cigar.data[j].n):
                        batch_info.sample_base_quals[m2] = batch_cigar.sample_base_quals_cigar.data[j].b
                        m2 += 1
                else:
                    m2 += batch_cigar.sample_base_quals_cigar.data[j].n

            # set ``mapqs_cigar``
            for j in range(batch_cigar.mapqs_cigar.size):

                if batch_cigar.mapqs_cigar.data[j].b != 0:
                    for _ in range(batch_cigar.mapqs_cigar.data[j].n):
                        batch_info.mapqs[m3] = batch_cigar.mapqs_cigar.data[j].b
                        m3 += 1
                else:
                    m3 += batch_cigar.mapqs_cigar.data[j].n

            # set ``read_pos_rank_cigar``
            for j in range(batch_cigar.read_pos_rank_cigar.size):

                if batch_cigar.read_pos_rank_cigar.data[j].b != 0:
                    for _ in range(batch_cigar.read_pos_rank_cigar.data[j].n):
                        batch_info.read_pos_rank[m4] = batch_cigar.read_pos_rank_cigar.data[j].b
                        m4 += 1
                else:
                    m4 += batch_cigar.read_pos_rank_cigar.data[j].n

            # set ``strands``
            for j in range(batch_cigar.strands_cigar.size):

                if batch_cigar.strands_cigar.data[j].b != '.':
                    for _ in range(batch_cigar.strands_cigar.data[j].n):
                        batch_info.strands[m5] = batch_cigar.strands_cigar.data[j].b
                        m5 += 1
                else:
                    m5 += batch_cigar.strands_cigar.data[j].n

        if self.position % 100000 == 0:
            logger.debug("Position %s:%s has %d base array, %d qual array, %d mapqs array, "
                         "%d pos-rank array, %d strands array." % (
                             self.chrid, self.position, total_base_array_size, total_qual_array_size,
                             total_mapqs_array_size, total_pos_rank_array_size, total_strand_array_size))

        return batch_info

    cdef _CigarStringArray _compress_string(self, char ** data, int size):
        cdef _CigarString *cigar = <_CigarString*> (calloc(size, sizeof(_CigarString)))

        cdef int last_index = 0
        cdef int i = 0, base_size
        for i in range(size):

            if i > 0:
                if strcmp(data[i], cigar[last_index].b) == 0:
                    # the two string is equal
                    cigar[last_index].n += 1
                else:
                    last_index += 1

                    base_size = strlen(data[i])
                    cigar[last_index].b = <char*> (calloc(base_size, sizeof(char)))
                    strcpy(cigar[last_index].b, data[i])

                    cigar[last_index].n = 1

            # i == 0
            else:
                last_index = 0

                base_size = strlen(data[i])
                cigar[last_index].b = <char*> (calloc(base_size, sizeof(char)))
                strcpy(cigar[last_index].b, data[i])
                cigar[last_index].n = 1

        cdef _CigarStringArray cigar_string_array
        cigar_string_array.size = last_index + 1
        cigar_string_array.data = <_CigarString*> (calloc(cigar_string_array.size, sizeof(_CigarString)))
        for i in range(cigar_string_array.size):
            cigar_string_array.data[i] = cigar[i]

        free(cigar)
        return cigar_string_array

    cdef _CigarCharArray _compress_char(self, char *data, int size):
        cdef _CigarChar *cigar = <_CigarChar*> (calloc(size, sizeof(_CigarChar)))

        cdef int last_index = 0
        cdef int i = 0
        for i in range(size):

            if i > 0:
                if data[i] == cigar[last_index].b:
                    cigar[last_index].n += 1
                else:
                    last_index += 1
                    cigar[last_index].b = data[i]
                    cigar[last_index].n = 1

            # i == 0
            else:

                last_index = 0

                cigar[last_index].b = data[i]
                cigar[last_index].n = 1

        cdef _CigarCharArray cigar_char_array
        cigar_char_array.size = last_index + 1
        cigar_char_array.data = <_CigarChar*> (calloc(cigar_char_array.size, sizeof(_CigarChar)))
        for i in range(cigar_char_array.size):
            cigar_char_array.data[i] = cigar[i]

        free(cigar)
        return cigar_char_array

    cdef _CigarIntArray _compress_int(self, int *data, int size):

        cdef _CigarInt *cigar = <_CigarInt*> (calloc(size, sizeof(_CigarInt)))

        cdef int last_index = 0
        cdef int i = 0
        for i in range(size):

            if i > 0:

                if data[i] == cigar[last_index].b:
                    cigar[last_index].n += 1
                else:
                    last_index += 1
                    cigar[last_index].b = data[i]
                    cigar[last_index].n = 1
            # i == 0
            else:

                last_index = 0
                cigar[last_index].b = data[i]
                cigar[last_index].n = 1

        cdef _CigarIntArray cigar_int_array

        cigar_int_array.size = last_index + 1
        cigar_int_array.data = <_CigarInt*> (calloc(cigar_int_array.size, sizeof(_CigarInt)))
        for i in range(cigar_int_array.size):
            cigar_int_array.data[i] = cigar[i]

        free(cigar)
        return cigar_int_array

cdef class BatchGenerator(object):
    """
    A class to generate batch informattion from a bunch of reads.
    """
    def __cinit__(self, bytes ref_name, long int reg_start, long int reg_end, FastaFile ref_fa, int sample_size,
                  options):
        """
        ``refresh``: mean the data in ``batch_heap`` will be changed.

        Constructor. Create a storage place for batchfile, and store the values of some flags which
        are used in the pysam CIGAR information.
        """
        # CIGAR here is the Mapping information in bwa.
        self.CIGAR_M = 0  # Match
        self.CIGAR_I = 1  # Insertion
        self.CIGAR_D = 2  # Deletion
        self.CIGAR_N = 3  # Skipped region from reference
        self.CIGAR_S = 4  # Soft clipping. Sequence is present in read
        self.CIGAR_H = 5  # Hard clipping. Sequence is not present in read
        self.CIGAR_P = 6  # Padding. Used for padded alignment
        self.CIGAR_EQ = 7  # Alignment match; sequence match
        self.CIGAR_X = 8  # Alignment match; sequence mismatch

        self.sample_size = sample_size
        self.ref_fa = ref_fa
        self.ref_name = ref_name
        self.reg_start = reg_start  # 1-base
        self.reg_end = reg_end  # 1-base

        self.ref_seq_start = max(0, self.reg_start - 200)
        self.ref_seq_end = min(self.reg_end + 200, self.ref_fa.references[self.ref_name].seq_length - 1)

        # initialization the BatchInfo for each position in `ref_name:reg_start-reg_end`
        cdef long int _pos  # `_pos` is 1-base system in the follow code.
        self.batch_heap = [BatchInfo(ref_name, _pos, self.ref_fa.get_character(self.ref_name, _pos-1), sample_size)
                           for _pos in range(reg_start, reg_end+1)]
        self.start_pos_in_batch_heap = reg_start  # 1-base, represent the first element in `batch_heap`
        self.options = options

        # the same definition with class ``BamReadBuffer`` in read.pyx
        self.filtered_read_counts_by_type = <int*> (calloc(7, sizeof(int)))
        if options.filter_duplicates == 0:
            self.filtered_read_counts_by_type[DUPLICATE] = -1

        if options.filter_reads_with_unmapped_mates == 0:
            self.filtered_read_counts_by_type[MATE_UNMAPPED] = -1

        if options.filter_reads_with_distant_mates == 0:
            self.filtered_read_counts_by_type[MATE_DISTANT] = -1

        if options.filterReadPairsWithSmallInserts == 0:
            self.filtered_read_counts_by_type[SMALL_INSERT] = -1

    def __dealloc__(self):
        """Clean up memory.
        """
        if self.filtered_read_counts_by_type != NULL:
            free(self.filtered_read_counts_by_type)

    cdef void create_batch_in_region(self, tuple region, cAlignedRead ** read_start, cAlignedRead ** read_end,
                                     int sample_index):  # The column index for the array in `batch_heap` represent a sample
        """Fetch batch information in a specific region."""
        cdef bytes chrom = region[0]
        cdef long int start = region[1]  # 1-base coordinate system
        cdef long int end = region[2]  # 1-base coordinate system

        if chrom != self.ref_name:
            logger.error("Error match chromosome (%s != %s)" % (chrom, self.ref_name))
            sys.exit(1)

        if self.start_pos_in_batch_heap < start or self.start_pos_in_batch_heap > end:
            logger.error("%s is not in region %s:%s-%s" % (
                self.start_pos_in_batch_heap, chrom, start, end))
            sys.exit(1)

        if start < self.ref_seq_start:
            logger.error("Start position (%s) is outside the reference region (%s)" % (
                start, "%s:%s-%s" % (self.ref_name, self.ref_seq_start, self.ref_seq_end)))
            sys.exit(1)

        if end > self.ref_seq_end:
            logger.error("End position (%s) is outside the reference region (%s)" % (
                end, "%s:%s-%s" % (self.ref_name, self.ref_seq_start, self.ref_seq_end)))
            sys.exit(1)

        cdef int read_num = 0
        while read_start != read_end:

            if Read_IsQCFail(read_start[0]):
                read_start += 1  # QC fail read move to the next one
                continue

            # still behind the region, do nothing but continue
            if read_start[0].end < start:
                read_start += 1
                continue

            # Break the loop when mapping start position is outside the region.
            if read_start[0].pos > end:
                break

            # get batch information here!
            self.get_batch_from_single_read_in_region(read_start[0], start, end, sample_index)

            read_num += 1  # how many reads in this regions
            read_start += 1  # move to the next read

        return

    cdef void get_batch_from_single_read_in_region(self, cAlignedRead *read, long int start, long int end,
                                                   int sample_index):  # The column index for the array in `batch_heap` which represent a sample
        """Check a single read for batch. 
        Batch positions are flagged by the CIGAR string. Pysam reports the CIGAR string information as a list 
        of tuples, where each tuple is a pair, and the first element gives the type of feature (match, insertion 
        or deletion), and the second element gives the number of nucleotides associated. 
        For example, [(0, 1), (1, 2), (0, 1)] is a 1 base match, a 2 base insertion, and a 1 base match.
        
        ``start``: 1-base
        ``end``: 1-base
        """
        cdef long int read_start_pos = read.pos  # 0-base

        cdef int ref_offset = 0
        cdef int read_offset = 0
        cdef int mapq = read.mapq
        cdef int cigar_length = read.cigar_len
        cdef char *read_seq = read.seq
        cdef char *read_qual = read.qual
        cdef bytes insert_seq = None
        cdef bytes deleted_seq = None
        cdef char map_strand = '-' if Read_IsReverse(read) else '+'

        cdef int cigar_flag = 0
        cdef int cigar_index = 0
        cdef int length = 0
        cdef long int ref_pos = 0
        cdef int pos_index, seq_index, base_index

        cdef BatchInfo batch_info
        for cigar_index in range(cigar_length):
            cigar_flag = read.cigar_ops[2 * cigar_index]
            length = read.cigar_ops[(2 * cigar_index) + 1]

            # An insertion take us further along the read, but not the reference
            if cigar_flag == self.CIGAR_I:

                if cigar_index > 0 and read.cigar_ops[(2 * cigar_index) - 2] == self.CIGAR_M:
                    pass
                elif cigar_index < cigar_length - 1 and read.cigar_ops[(2 * cigar_index) + 2] == self.CIGAR_M:
                    pass
                else:
                    read_offset += length
                    continue

                insert_seq = read_seq[read_offset:read_offset + length]
                ref_pos = read_start_pos + ref_offset - 1  # `ref_pos` is 0-base
                ref_pos += 1  # `ref_pos` is 1-base

                # do not use insertion with Ns in them
                if (start <= ref_pos <= end) and (insert_seq.count("N") == 0):

                    pos_index = ref_pos - self.start_pos_in_batch_heap
                    batch_info = self.batch_heap[pos_index]

                    if batch_info.is_empty[sample_index]:
                        # Just record the information of first read, so do not update if it's not empty
                        batch_info.update_info_by_index(sample_index, self.ref_name, ref_pos, mapq, map_strand,
                                                        "+%s" % insert_seq, 0, read_offset)

                read_offset += length

            # A deletion take us further along the reference, but not the read
            elif cigar_flag == self.CIGAR_D:
                if cigar_index > 0 and read.cigar_ops[(2 * cigar_index) - 2] == self.CIGAR_M:
                    pass
                elif cigar_index < cigar_length - 1 and read.cigar_ops[(2 * cigar_index) + 2] == self.CIGAR_M:
                    pass
                else:
                    ref_offset += length
                    continue

                deleted_seq = self.ref_fa.get_sequence(self.ref_name, read_start_pos + ref_offset,
                                                       read_start_pos + ref_offset + length)
                ref_pos = read_start_pos + ref_offset - 1
                ref_pos += 1  # `ref_pos` is 1-base
                # do not use deletion with Ns in them
                if (start <= ref_pos <= end) and (deleted_seq.upper().count("N") == 0):

                    pos_index = ref_pos - self.start_pos_in_batch_heap
                    batch_info = self.batch_heap[pos_index]

                    if batch_info.is_empty[sample_index]:
                        # Just record the information of first read, so do not update if it's not empty
                        batch_info.update_info_by_index(sample_index, self.ref_name, ref_pos, mapq, map_strand,
                                                        "-%s" % deleted_seq, 0, read_offset)

                ref_offset += length

            # A match take us further along the reference and the read
            elif cigar_flag == self.CIGAR_M or cigar_flag == self.CIGAR_EQ or cigar_flag == self.CIGAR_X:

                self._get_matchbase_from_read_segment(
                    sample_index, read_seq, read_qual, mapq, map_strand, read_start_pos,
                    read_offset, ref_offset, length, start, end
                )

                ref_offset += length
                read_offset += length

            elif cigar_flag == self.CIGAR_N:
                ref_offset += length

            elif cigar_flag == self.CIGAR_S:
                read_offset += length
                # We have to move back read position when there is a soft-clipping at the
                # beginning of reads
                if cigar_index == 0:
                    ref_offset += length

            # Hard clipping. Sequence is not present in read.
            elif cigar_flag == self.CIGAR_H:
                continue

            # Padding. We do nothing here.
            elif cigar_flag == self.CIGAR_P:
                continue

            # Other kinds of flag. Just don't care about them.
            else:
                continue

        return

    cdef void _get_matchbase_from_read_segment(self, int sample_index,
                                               char *read_seq,
                                               char *read_qual,
                                               int mapq,
                                               char map_strand,
                                               int read_start,
                                               int read_offset,
                                               int ref_offset,
                                               int seglength,
                                               long int region_start,
                                               long int region_end):
        """Get BatchElement from a particular CIGAR_M segment in read.
        
        Parameter
        =========
            ``read_seq``: The complete sequence of bases in the read
            ``read_qual``: The complete sequence of quality scores in the read
            ``read_start``: The starting position, in the reference sequence, of the read
            ``read_offset``: int
                If we're not starting from the beginning of the read, then how far along 
                the read sequence to start.
            ``ref_offset``: int
                If the read contains indels, we need to offset our reference position accordingly, 
                by this much.
        """
        # ignore all the positions which are not in [region_start, region_end]
        if (read_start + ref_offset > region_end) or (read_start + ref_offset + seglength < region_start):
            return

        cdef int base_qual = 0
        cdef int base_index = 0
        cdef int pos_index = 0
        cdef char base_char[2]

        cdef int index = 0
        cdef long int ref_pos = 0
        cdef BatchInfo batch_info
        for index in range(seglength):

            ref_pos = read_start + ref_offset + index
            ref_pos += 1  # `ref_pos` is 1-base

            if ref_pos < region_start:
                continue

            if ref_pos > region_end:
                break

            base_index = read_offset + index
            base_char[0] = read_seq[base_index]
            base_char[1] = '\0'
            base_qual = read_qual[base_index]

            pos_index = ref_pos - self.start_pos_in_batch_heap
            batch_info = self.batch_heap[pos_index]

            if batch_info.is_empty[sample_index]:
                # Just record the information of first read, so do not update if it's not empty
                batch_info.update_info_by_index(sample_index, self.ref_name, ref_pos, mapq, map_strand,
                                                base_char, base_qual, base_index)

        return
