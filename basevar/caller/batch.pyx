# cython: profile=True
"""
Author: Shujia Huang
Date: 2019-06-05 10:59:21
"""
from libc.stdio cimport fprintf, stderr, stdout
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdlib cimport calloc, realloc, free
from libc.string cimport strcpy, strcmp, strlen

from basevar.utils cimport c_min, c_max
from basevar.io.read cimport cAlignedRead
from basevar.io.htslibWrapper cimport Read_IsQCFail
from basevar.io.htslibWrapper cimport Read_IsReverse
from basevar.datatype.genomeregion cimport GenomeRegion
from basevar.datatype.dynamicstring cimport dstring, dstring_init, dstring_append, \
    dstring_append_char, dstring_append_long, dstring_append_int

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
    def __cinit__(self, char *chrid, long int position, char *ref_base, int size):
        self.size = size  # default size is as the same as capacity
        self.__capacity = size

        self.chrid = <char *>calloc(strlen(chrid)+1, sizeof(char))
        strcpy(self.chrid, chrid)

        self.position = position

        self.ref_base = <char *>calloc(strlen(ref_base)+1, sizeof(char))
        strcpy(self.ref_base, ref_base)

        self.depth = 0
        self.strands = <char*> (calloc(self.__capacity, sizeof(char)))
        self.sample_bases = <char**> (calloc(self.__capacity, sizeof(char*)))
        self.sample_base_quals = <int*> (calloc(self.__capacity, sizeof(int)))
        self.mapqs = <int*> (calloc(self.__capacity, sizeof(int)))
        self.read_pos_rank = <int*> (calloc(self.__capacity, sizeof(int)))
        self.is_empty = <int*> (calloc(self.__capacity, sizeof(int)))

        if (self.strands == NULL) or (self.sample_bases == NULL) or (self.sample_base_quals == NULL) \
                or (self.mapqs == NULL) or (self.read_pos_rank == NULL) or (self.is_empty == NULL):
            fprintf(stderr, "[ERROR] Could not allocate memory for self.samples_info in BatchInfo.\n")
            exit(EXIT_FAILURE)

        self.set_empty()

    cdef void set_size(self, size_t size):
        """call this function to adjust the self.size, if the real size is not equal to capacity"""
        if size > self.__capacity:
            fprintf(stderr, "[ERROR] The size (%zu) is larger than the capacity(%zu) in ``BatchInfo``.", size,
                    self.__capacity)
            exit(EXIT_FAILURE)

        self.size = size
        return

    cdef void set_empty(self):
        # initialization and set empty mark for all the element, 1=>empty, 0=> not empty
        cdef int i
        for i in range(self.__capacity):
            self.is_empty[i] = 1
            self.mapqs[i] = 0
            self.strands[i] = '.'

            self.sample_bases[i] = <char*>calloc(2, sizeof(char))
            strcpy(self.sample_bases[i], "N")

            self.sample_base_quals[i] = 0
            self.read_pos_rank[i] = 0

        self.depth = 0
        return

    cdef dstring get_str(self):

        cdef dstring ds
        dstring_init(&ds, 128)
        fprintf(stdout, ">>>>>>>>>>>> %d, %d, %s\n", ds.size, ds.__capacity, ds.s)

        dstring_append(&ds, self.chrid)
        dstring_append_char(&ds, '\t')
        dstring_append_long(&ds, self.position)
        fprintf(stdout, "*******: %d, %d, %s\n", ds.size, ds.__capacity, ds.s)
        dstring_append_char(&ds, '\t')
        dstring_append(&ds, self.ref_base)
        dstring_append_char(&ds, '\t')
        dstring_append_long(&ds, self.depth)
        dstring_append_char(&ds, '\t')

        cdef unsigned i = 0
        if self.depth > 0:

            for i in range(self.size):
                if i > 0:
                    dstring_append_char(&ds, ',')
                dstring_append_int(&ds, self.mapqs[i])
            dstring_append_char(&ds, '\t')

            for i in range(self.size):
                if i > 0:
                    dstring_append_char(&ds, ',')
                dstring_append(&ds, self.sample_bases[i]) # Must all be upper
            dstring_append_char(&ds, '\t')

            for i in range(self.size):
                if i > 0:
                    dstring_append_char(&ds, ',')
                dstring_append_int(&ds, self.sample_base_quals[i])
            dstring_append_char(&ds, '\t')

            for i in range(self.size):
                if i > 0:
                    dstring_append_char(&ds, ',')
                dstring_append_int(&ds, self.read_pos_rank[i])
            dstring_append_char(&ds, '\t')

            for i in range(self.size):
                if i > 0:
                    dstring_append_char(&ds, ',')
                dstring_append_char(&ds, self.strands[i])
        else:
            dstring_append(&ds, ".\t.\t.\t.\t.")

        return ds

    cdef void update_info_by_index(self, int index, char *_target_chrom, unsigned long _target_position, int mapq,
                                   char map_strand, char *read_base, int base_qual, int read_pos_rank):
        """Update information"""
        # fprintf (stdout, "---- before update: %s, %ld, %d, %s\n", self.chrid, self.position, self.depth, read_base)
        if strcmp(self.chrid, _target_chrom) != 0 or (_target_position != self.position):
            fprintf(stderr, "[ERROR] Chromosome(%s, %s) or position(%ld, %ld) not exactly match.\n",
                self.chrid, _target_chrom, self.position, _target_position)

        if read_base[0] != '-' and read_base[0] != '+':
            self.depth += 1


        if self.is_empty[index]:
            self.is_empty[index] = 0  # not empty

        self.mapqs[index] = mapq
        self.strands[index] = map_strand
        self.sample_base_quals[index] = base_qual
        self.read_pos_rank[index] = read_pos_rank

        free(self.sample_bases[index])
        self.sample_bases[index] = <char*> (calloc(strlen(read_base)+1, sizeof(char)))
        strcpy(self.sample_bases[index], read_base)

        fprintf (stdout, "---- After update: %s, %ld, %d, %s, %d, %s\n", self.chrid, self.position, self.depth, read_base, strlen(read_base), self.sample_bases[index])
        return

    cdef void clear(self):
        self.__dealloc__()
        return

    def __dealloc__(self):
        """Free memory"""
        free(self.chrid)
        free(self.ref_base)

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
    def __cinit__(self, char *chrid, long int position, char *ref_base, unsigned long array_size):
        """Allocate an array of size 'size', with initial values 'init'.
        """
        # base information at specific position!
        self.chrid = chrid
        self.position = position
        self.ref_base = ref_base

        self.array = <BatchCigar*>calloc(array_size, sizeof(BatchCigar))
        if self.array == NULL:
            fprintf(stderr, "[ERROR] Could not allocate PositionBatchCigarArray!!\n")
            exit(EXIT_FAILURE)

        self.__size = 0  # We don't put anything in here yet
        self.__capacity = array_size
        self.__sample_number = 0  # The number of samples which have been store in this array
        self.__depth = 0

    def __dealloc__(self):
        """
        Free memory
        """
        cdef BatchCigar batch_cigar
        cdef unsigned i = 0
        if self.array != NULL:
            for i in range(self.__size):
                batch_cigar = self.array[i]
                free(batch_cigar.sample_bases_cigar.data)
                free(batch_cigar.sample_base_quals_cigar.data)
                free(batch_cigar.read_pos_rank_cigar.data)
                free(batch_cigar.mapqs_cigar.data)
                free(batch_cigar.strands_cigar.data)

            free(self.array)

    cdef unsigned long size(self):
        return self.__size

    cdef void append(self, BatchInfo value):

        if value.chrid != self.chrid or value.position != self.position:
            fprintf(stderr,
                    "[ERROR] Error in ``PositionBatchCigarArray`` when call append()! "
                    "Chromosome(%s, %s) or position(%ld, %ld) not exactly match \n",
                    self.chrid, value.chrid, self.position, value.position)
            exit(EXIT_FAILURE)

        cdef BatchCigar *temp = NULL
        if self.__size == self.__capacity:
            temp = <BatchCigar*> (realloc(self.array, 2 * sizeof(BatchCigar) * self.__capacity))

            if temp == NULL:
                fprintf(stderr, "[ERROR] Could not re-allocate PositionBatchCigarArray!!\n")
                exit(EXIT_FAILURE)
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
            fprintf(stdout,
                    "[INFO] Position %s:%ld has %d base array, %d qual array, %d mapqs array, "
                    "%d pos-rank array, %d strands array.\n", self.chrid, self.position, total_base_array_size,
                    total_qual_array_size, total_mapqs_array_size, total_pos_rank_array_size, total_strand_array_size)

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
    A class to generate batch information from a bunch of reads.
    """
    def __cinit__(self, const GenomeRegion region,  # The coordinate system of ``region`` is 1-base
                  FastaFile ref_fa,
                  int sample_size,
                  BaseTypeCmdOptions options):
        """
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
        self.CIGAR_EQ = 7 # Alignment match; sequence match
        self.CIGAR_X = 8  # Alignment match; sequence mismatch

        self.region = region  # 1-base system
        self.sample_size = sample_size
        self.ref_fa = ref_fa

        self.ref_seq_start = c_max(0, self.region.start-200)
        self.ref_seq_end = c_min(self.region.end+200, self.ref_fa.references[self.region.chrom].seq_length-1)

        # initialization the BatchInfo for each position in `ref_name:reg_start-reg_end`
        cdef unsigned long _pos  # `_pos` is 1-base system in the follow code.
        self.batch_heap = [BatchInfo(self.region.chrom, _pos, self.ref_fa.get_character(self.region.chrom, _pos-1), sample_size)
                           for _pos in range(self.region.start, self.region.end+1)]
        self.start_pos_in_batch_heap = self.region.start  # 1-base, represent the first element in `batch_heap`
        self.options = options

        # the same definition with class ``BamReadBuffer`` in read.pyx
        self.filtered_read_counts_by_type = <int *>calloc(7, sizeof(int))
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

    cdef void create_batch_in_region(self, const GenomeRegion region, cAlignedRead **read_start, cAlignedRead **read_end,
                                     int sample_index):  # The column index for the array in `batch_heap` represent a sample
        """Fetch batch information in a specific region."""
        if strcmp(self.region.chrom, region.chrom) != 0:
            fprintf(stderr, "[ERROR] Error match chromosome (%s != %s)\n", self.region.chrom, region.chrom)
            exit(EXIT_FAILURE)

        if self.start_pos_in_batch_heap < region.start or self.start_pos_in_batch_heap > region.end:
            fprintf(stderr, "[ERROR] %ld is not in region %s:%ld-%ld\n",
                self.start_pos_in_batch_heap, region.chrom, region.start, region.end)
            exit(EXIT_FAILURE)

        if region.start < self.ref_seq_start:
            fprintf(stderr, "[ERROR] Start position (%ld) is outside the reference region (%s:%ld-%ld)\n",
                region.start, self.region.chrom, self.ref_seq_start, self.ref_seq_end)
            exit(EXIT_FAILURE)

        if region.end > self.ref_seq_end:
            fprintf(stderr, "[ERROR] End position (%ld) is outside the reference region (%s:%ld-%ld)\n",
                region.end, self.region.chrom, self.ref_seq_start, self.ref_seq_end)
            exit(EXIT_FAILURE)

        cdef int read_num = 0
        cdef long int read_start_pos  # delete
        while read_start != read_end:

            if Read_IsQCFail(read_start[0]):
                read_start += 1  # QC fail read move to the next one
                continue

            # still behind the region, do nothing but continue
            if read_start[0].end < region.start:
                read_start += 1
                continue

            # Break the loop when mapping start position is outside the region.
            if read_start[0].pos > region.end:
                break

            # get batch information here!
            self.get_batch_from_single_read_in_region(read_start[0], region.start, region.end, sample_index)
            read_start_pos = read_start[0].pos
            read_num += 1  # how many reads in this regions
            fprintf(stdout, "Read:[%ld, %ld], %ld,  read_num: %d\n", region.start, region.end, read_start_pos, read_num)
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
                        batch_info.update_info_by_index(sample_index, self.region.chrom, ref_pos, mapq, map_strand,
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

                deleted_seq = self.ref_fa.get_sequence(self.region.chrom, read_start_pos + ref_offset,
                                                       read_start_pos + ref_offset + length)
                ref_pos = read_start_pos + ref_offset - 1
                ref_pos += 1  # `ref_pos` is 1-base
                # do not use deletion with Ns in them
                if (start <= ref_pos <= end) and (deleted_seq.upper().count("N") == 0):

                    pos_index = ref_pos - self.start_pos_in_batch_heap
                    batch_info = self.batch_heap[pos_index]

                    if batch_info.is_empty[sample_index]:
                        # Just record the information of first read, so do not update if it's not empty
                        batch_info.update_info_by_index(sample_index, self.region.chrom, ref_pos, mapq, map_strand,
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
                                               const char *read_seq,
                                               const char *read_qual,
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

            # if batch_info.is_empty[sample_index]:
            if batch_info.is_empty[sample_index]:
                # Just record the information of first read, so do not update if it's not empty
                batch_info.update_info_by_index(sample_index, self.region.chrom, ref_pos, mapq, map_strand,
                                                base_char, base_qual, base_index)

        return
