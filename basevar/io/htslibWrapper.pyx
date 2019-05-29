"""Wrapper for htslib"""
import os

from basevar.log import logger

# valid types for sam headers
VALID_HEADER_TYPES = {"HD": dict,
                      "SQ": list,
                      "RG": list,
                      "PG": list,
                      "CO": list}

# order of records within sam headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO")

# type conversions within sam header records
VALID_HEADER_FIELDS = {"HD": {"VN": str, "SO": str, "GO": str},
                       "SQ": {"SN": str, "LN": int, "AH": str, "AS": str, "M5": str, "UR": str, "SP": str},
                       "RG": {"ID": str, "SM": str, "LB": str, "DS": str, "PU": str, "PI": str, "CN": str, "DT": str,
                              "PL": str, "PG": str},
                       "PG": {"ID": str, "VN": str, "CL": str, "PN": str, "PP": str, "DS": str}, }

# output order of fields within records
VALID_HEADER_ORDER = {"HD": ("VN", "SO", "GO"),
                      "SQ": ("SN", "LN", "AH", "AS", "M5", "UR", "SP"),
                      "RG": ("ID", "SM", "LB", "DS", "PU", "PI", "CN", "DT", "PL", "PG"),
                      "PG": ("ID", "VN", "CL"), }

######################################################################
## Public methods
######################################################################
cdef class Samfile:
    """The class for SAM/BAM/CRAM file.

    *(filename, mode='r', referencenames = None, referencelengths = None, text = NULL, header = None)*

    A *SAM* file. The file is automatically opened.

    *mode* should be ``r`` for reading or ``w`` for writing. The default is text mode so for binary
    (:term:`BAM`) I/O you should append ``b`` for compressed or ``u`` for uncompressed :term:`BAM` output.
    Use ``h`` to output header information  in text (:term:`TAM`)  mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``.
    Currently valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb`` and ``wbu``.

    so to open a :term:`BAM` file for reading::

        f=Samfile('ex1.bam','rb')


    For writing, the header of a :term:`TAM` file/:term:`BAM` file can be constituted from several
    sources:

        2. If *header* is given, the header is build from a multi-level dictionary.
        The first level are the four types ('HD', 'SQ', ...). The second level is then a list of lines,
        with each line being a list of tag-value pairs.

        3. If *text* is given, new header text is copied from raw text.

        4. The names (*referencenames*) and lengths (*referencelengths*) are supplied directly as lists.

    If an index for a BAM file exists (.bai), it will be opened automatically. Without an index random
    access to reads via :meth:`fetch` and :meth:`pileup` is disabled.
    """
    def __cinit__(self, filename):
        """Constructor.
        """
        self.filename = filename
        self.samfile = NULL
        self.the_header = NULL
        self.index = NULL
        self.lock = None

    def __dealloc__(self):
        """Clean up.
        Here I've pasted the code from several clean-up
        methods, as calling methods is not guaranteed to work in
        this function.
        """
        # remember: dealloc cannot call other methods
        # Note that __del__ is not called.
        if self.samfile != NULL:
            sam_close(self.samfile)
            self.samfile = NULL

        if self.index != NULL:
            hts_idx_destroy(self.index)
            self.index = NULL

        if self.the_header != NULL:
            bam_hdr_destroy(self.the_header)
            self.the_header = NULL

    cdef void clear_header(self):
        """ Clear all the header data. """
        if self.the_header != NULL:
            bam_hdr_destroy(self.the_header)
            self.the_header = NULL

    cdef void clear_index(self):
        """Clear all the index file data. """
        print("CLEARING INDEX!!!!!!")
        if self.index != NULL:
            hts_idx_destroy(self.index)
            self.index = NULL

    cdef int _is_bam(self):
        return self.samfile.is_bin

    cdef int _is_cram(self):
        return self.samfile.is_cram

    cdef int _is_open(self):
        """return true if samfile has been opened."""
        return self.samfile != NULL

    cdef bool _has_index(self):
        """return true if samfile has an existing (and opened) index."""
        return self.index != NULL

    cdef void _open(self, mode, load_index=True):
        """open a sam/bam/cram file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
        """
        if mode not in ("r", "rb", "rbh"):
            raise StandardError, "invalid file opening mode `%s`" % mode

        if self.samfile != NULL:
            if load_index and self.index == NULL:
                self.index = sam_index_load(self.samfile, self.filename)
                if self.index == NULL:
                    raise IOError("Error while opening index for file `%s`. "
                                  "Check that index exists " % self.filename)
            return

        if mode[0] == "r":
            self._open_bamfile(mode)
        else:
            raise StandardError, "BAM file is read-only."

        if self.samfile == NULL:
            raise IOError("Could not open file `%s`. Check that file/path exists." % self.filename)

        if self._is_bam() or self._is_cram():
            # returns NULL if there is no index or index could not be opened
            if load_index and self.index == NULL:
                self.index = sam_index_load(self.samfile, self.filename)
                if self.index == NULL:
                    raise IOError("Error while opening index for file `%s`. "
                                  "Check that index exists " % self.filename)

    cdef void _open_bamfile(self, mode):
        """Open BamFile
        """
        if os.path.exists(self.filename):
            self.samfile = sam_open(self.filename, mode)
            self.the_header = bam_hdr_init()
            self.the_header = sam_hdr_read(self.samfile)
        else:
            self.samfile = NULL
            self.the_header = NULL

        return

    cdef char* getrname(self, int tid):
        """Convert numerical :term:`tid` into :ref:`reference` name."""
        if not 0 <= tid < self.the_header.n_targets:
            raise ValueError("tid (%s) out of range 0<=tid<%i" % (tid, self.the_header.n_targets))

        return self.the_header.target_name[tid]

    cpdef ReadIterator fetch(self, const char *region):
        """Fetch reads from a specified region."""

        # Need to load the index for new queries.
        if not self._is_open():
            self._open('r', loadIndex=True)

        if self._is_bam() or self._is_cram():
            return ReadIterator(self, region)
        else:
            raise StandardError, "Random access query only allowed for BAM/CRAM files."

    cpdef void close(self):
        """closes file."""

        if not self._is_cram():
            if self.samfile != NULL:
                sam_close(self.samfile)
                self.samfile = NULL

            if self.index != NULL:
                hts_idx_destroy(self.index)
                self.index = NULL

            if self.the_header != NULL:
                bam_hdr_destroy(self.the_header)
                self.the_header = NULL

    property nreferences:
        """number of :term:`reference` sequences in the file."""

        def __get__(self):
            return self.the_header.n_targets

    property references:
        """tuple with the names of :term:`reference` sequences."""
        def __get__(self):
            t = []
            cdef int x = 0
            for x in range(self.the_header.n_targets):
                t.append(self.the_header.target_name[x])
            return tuple(t)

    property lengths:
        """
        tuple of the lengths of the :term:`reference` sequences. The lengths are 
        in the same order as :attr:`pysam.Samfile.reference`
        """
        def __get__(self):
            t = []
            cdef int x = 0
            for x in range(self.the_header.n_targets):
                t.append(self.the_header.target_len[x])
            return tuple(t)

    property text:
        """full contents of the :term:`sam file` header as a string."""

        def __get__(self):
            # create a temporary 0-terminated copy
            cdef char *t
            t = <char*> calloc(self.the_header.l_text + 1, sizeof(char))
            memcpy(t, self.the_header.text, self.the_header.l_text)
            cdef bytes result = t

            free(t)
            return result

    property header:
        """header information within the :term:`sam file`. 
        The records and fields are returned as a two-level dictionary.
        """
        def __get__(self):
            result = {}

            if self.the_header.text != NULL:
                # convert to python string (note: call self.text to create 0-terminated string)
                t = self.text
                for line in t.split("\n"):
                    if not line.strip(): continue

                    if not line.startswith("@"):
                        raise StandardError, "Header line without '@': '%s. Total header text is %s'" % (line, t)

                    fields = line[1:].split("\t")
                    record = fields[0]
                    assert record in VALID_HEADER_TYPES, "header line with invalid type '%s': '%s'" % (record, line)

                    # treat comments
                    if record == "CO":
                        if record not in result: result[record] = []
                        result[record].append("\t".join(fields[1:]))
                        continue

                    # the following is clumsy as generators do not work?
                    x = {}
                    for field in fields[1:]:
                        key, value = field.split(":", 1)
                        if key not in VALID_HEADER_FIELDS[record]:
                            raise ValueError("unknown field code '%s' in record '%s'" % (key, record))

                        x[key] = VALID_HEADER_FIELDS[record][key](value)

                    if VALID_HEADER_TYPES[record] == dict:
                        if record in result:
                            raise ValueError("multiple '%s' lines are not permitted" % record)
                        result[record] = x

                    elif VALID_HEADER_TYPES[record] == list:
                        if record not in result:
                            result[record] = []

                        result[record].append(x)

            return result

cdef class ReadIterator:
    """
    Iterates over mapped reads in a region.
    """

    def __cinit__(self, Samfile samfile, const char *region):
        """
        Constructor.
        """

        self.theSamfile = NULL
        self.theIterator = NULL
        self.b = NULL

        if not samfile._is_open():
            raise StandardError, "Samfile %s is not open. Cannot read from file." % (samfile.filename)

        self.theSamfile = samfile.samfile

        # Load from BAM by querying index for file off-set
        if samfile._has_index():
            self.theIterator = sam_itr_querys(samfile.index, samfile.the_header, region)
        else:
            raise StandardError, ("Cannot retrieve random region from Samfile %s, as it "
                                  "does not have an index" % samfile.filename)

        self.b = bam_init1()

    def __dealloc__(self):
        """remember: dealloc cannot call other methods!"""

        if self.the_iterator != NULL:
            sam_itr_destroy(self.the_iterator)

        if self.b != NULL:
            bam_destroy1(self.b)

    cdef cAlignedRead* get(self, int store_rgID, char** rgID):
        cdef bam1_core_t* c = &self.b.core
        cdef uint8_t* s = bam_get_seq(self.b)
        cdef uint8_t* q = bam_get_qual(self.b)
        cdef int len_seq = c.l_qseq

        if len_seq == 0:
            return NULL

        if q[0] == 0xff:
            return NULL

        cdef cAlignedRead* the_read = <cAlignedRead*> malloc(sizeof(cAlignedRead))
        cdef char* seq = <char*> malloc((len_seq + 1) * sizeof(char))
        cdef char* qual = <char*> malloc((len_seq + 1) * sizeof(char))

        assert the_read != NULL
        assert seq != NULL
        assert qual != NULL

        # Try to grab the read-group tag value
        cdef uint8_t* v = NULL
        cdef char* temp_rgID = NULL
        cdef int len_rgID = 0

        if store_rgID:
            v = bam_aux_get(self.b, "RG")
            if v != NULL:
                temp_rgID = bam_aux2Z(v)
                len_rgID = strlen(temp_rgID)

                rgID[0] = <char*> (calloc(len_rgID + 1, sizeof(char)))
                strcpy(rgID[0], temp_rgID)
            else:
                rgID[0] = NULL

        cdef int i = 0
        for i in range(len_seq):
            seq[i] = self._getBase(s, i)
            qual[i] = q[i]
            assert qual[i] <= 93
            assert qual[i] >= 0

        seq[lenSeq] = '\0'
        qual[lenSeq] = '\0'

        readStart = c.pos
        cdef short*cigarOps = <short*> malloc(2 * c.n_cigar * sizeof(short))
        assert cigarOps != NULL
        cdef uint32_t *cigar = bam_get_cigar(self.b)
        for i from 0 <= i < c.n_cigar:
            cigarFlag = bam_cigar_op(cigar[i])
            cigarFlagLen = bam_cigar_oplen(cigar[i])
            cigarOps[2 * i] = cigarFlag
            cigarOps[(2 * i) + 1] = cigarFlagLen

            # Soft-clipping of sequence at start of read changes the mapping
            # position. Recorded mapping pos is that of the first aligned (not soft-clipped)
            # base. I want to adjust this so that the read start refers to the first base in
            # the read.
            if i == 0 and cigarFlag == 4:
                readStart -= cigarFlagLen

        theRead.seq = seq
        theRead.qual = qual
        theRead.cigarOps = cigarOps
        theRead.hash = NULL
        theRead.mateChromID = c.mtid
        theRead.cigarLen = c.n_cigar
        theRead.chromID = c.tid
        theRead.rlen = lenSeq
        theRead.pos = readStart
        theRead.end = bam_endpos(self.b)
        theRead.insertSize = c.isize
        theRead.matePos = c.mpos
        theRead.bitFlag = c.flag
        theRead.mapq = c.qual

        Read_SetUnCompressed(theRead)

        return theRead

    cdef int cnext(self) nogil:
        '''
        cversion of iterator. Used by IteratorColumn
        '''
        return sam_itr_next(self.theSamfile, self.theIterator, self.b) >= 0

    cdef char _getBase(self, uint8_t *s, int i):
        cdef char* baseLookup = "=ACMGRSVTWYHKDBN"
        return baseLookup[bam_seqi(s, i)]


cdef void destroy_read(cAlignedRead*theRead):
    """
    De-allocate memory for read.
    """

    if theRead.seq != NULL:
        free(theRead.seq)

    if theRead.qual != NULL:
        free(theRead.qual)

    if theRead.cigarOps != NULL:
        free(theRead.cigarOps)

    if theRead.hash != NULL:
        free(theRead.hash)

    free(theRead)


cdef void compress_seq(cAlignedRead*read, char*refSeq):
    """
    Does exactly what it says on the tin.
    """
    cdef char*seq = read.seq
    cdef char*newSeq = <char*> (alloca(sizeof(char) * read.rlen * 2))
    cdef char*finalSeq = NULL
    cdef int i = 0
    cdef int nMatches = 0
    cdef int newSeqIndex = 0
    cdef int totalMatches = 0

    for i in range(read.rlen):

        # Ref match. Store count if > 40
        if seq[i] == refSeq[i]:
            totalMatches += 1

            if nMatches == 40:
                newSeq[newSeqIndex] = nMatches
                nMatches = 0
                newSeqIndex += 1

            nMatches += 1

        # Ref mis-match. Store base.
        else:
            if nMatches > 0:
                newSeq[newSeqIndex] = nMatches
                nMatches = 0
                newSeqIndex += 1

            newSeq[newSeqIndex] = seq[i]
            newSeqIndex += 1

    # If we finish off with a load of matches, catch those
    if nMatches > 0:
        newSeq[newSeqIndex] = nMatches
        nMatches = 0
        newSeqIndex += 1

    newSeq[newSeqIndex] = 0
    finalSeq = <char*> malloc(sizeof(char) * (newSeqIndex + 1))
    strcpy(finalSeq, newSeq)
    finalSeq[newSeqIndex] = 0
    free(read.seq)
    #free(newSeq)
    read.seq = finalSeq

###################################################################################################

cdef void compress_qual(cAlignedRead*read, int qualBinSize):
    """
    Does exactly what it says on the tin.
    """
    cdef int i = 0
    cdef char*qual = read.qual
    #cdef char* newQual = <char*>(malloc(sizeof(char)*2*read.rlen))
    cdef char*newQual = <char*> (alloca(sizeof(char) * 2 * read.rlen))
    cdef char*finalQual = NULL
    cdef char lastChar = 0
    cdef int lastCount = 0
    cdef int newQualIndex = 0

    if qualBinSize > 1:
        for i in range(read.rlen):
            qual[i] = (qual[i] / qualBinSize) * qualBinSize

    for i in range(read.rlen):

        assert qual[i] >= 0, "Shit 3. Qual = %s" % (qual[i])
        assert qual[i] <= 93, "Shit 4!. Qual = %s" % (qual[i])

        if i == 0:
            newQual[newQualIndex] = qual[i] + 33
            newQualIndex += 1
            lastChar = qual[i]
            lastCount = 1
        else:
            if qual[i] == lastChar:
                lastCount += 1
            else:
                newQual[newQualIndex] = lastCount
                newQualIndex += 1
                newQual[newQualIndex] = qual[i] + 33
                newQualIndex += 1
                lastChar = qual[i]
                lastCount = 1

    if lastCount > 0:
        newQual[newQualIndex] = lastCount
        newQualIndex += 1

    finalQual = <char*> malloc(sizeof(char) * (newQualIndex + 1))
    newQual[newQualIndex] = 0
    strcpy(finalQual, newQual)
    finalQual[newQualIndex] = 0
    #free(newQual)
    free(read.qual)
    read.qual = finalQual

cdef void uncompress_seq(cAlignedRead*read, char*refSeq):
    """
    Does exactly what it says on the tin.
    """
    cdef int refIndex = 0
    cdef int newSeqIndex = 0
    cdef int i = 0
    cdef int j = 0
    cdef int seqLen = strlen(read.seq)
    cdef char*seq = read.seq
    cdef char*newSeq = <char*> (malloc(sizeof(char*) * (read.rlen + 1)))

    for i in range(seqLen):
        if seq[i] <= 40:
            for j in range(seq[i]):
                newSeq[newSeqIndex] = refSeq[j + refIndex]
                newSeqIndex += 1
            refIndex += seq[i]
        else:
            newSeq[newSeqIndex] = seq[i]
            refIndex += 1
            newSeqIndex += 1

    free(read.seq)
    read.seq = newSeq
    read.seq[read.rlen] = 0

cdef void uncompress_qual(cAlignedRead*read):
    """
    Does exactly what it says on the tin.
    """
    cdef char*newQual = <char*> (malloc(sizeof(char) * (read.rlen + 1)))
    cdef char*qual = read.qual

    cdef int lenQual = strlen(read.qual)
    cdef int i = 0
    cdef int j = 0
    cdef int newQualIndex = 0

    #logger.info("Length = %s or %s" %(lenQual, len([x for x in qual])))
    #logger.info(",".join( [str(x) for x in qual] ) )

    for i in range(0, lenQual - 1, 2):
        for j in range(qual[i + 1]):
            newQual[newQualIndex] = qual[i] - 33
            assert qual[i] - 33 >= 0, "Shit 1. Qual = %s. Count = %s!" % (qual[i], qual[i + 1])
            assert qual[i] - 33 <= 93, "Shit 2!. Qual = %s. Count = %s" % (qual[i], qual[i + 1])
            newQualIndex += 1

    assert read.rlen == newQualIndex

    free(read.qual)
    read.qual = newQual
    read.qual[read.rlen] = 0

cdef void compress_read(cAlignedRead*read, char*refSeq, int refStart, int refEnd, int qualBinSize, int fullComp):
    """To save memory use, we compress reads. The sequence is compressed using reference-based compression.
    
    The qualities are compressed using zlib and (optionally) using course binning. We delete the read-group tag, 
    and the cigar string information (this can only be done after the read has been used for candidate generation).
    A bit is set in the bit-field to label the read as compressed or uncompressed.
    """
    if Read_IsCompressed(read):
        return

    if read.seq != NULL:
        compressSeq(read, refSeq + (read.pos - refStart))

    if read.qual != NULL:
        compressQual(read, qualBinSize)

    if read.hash != NULL:
        free(read.hash)
        read.hash = NULL

    Read_SetCompressed(read)

cdef void uncompress_read(cAlignedRead*read, char*refseq, int refstart, int refend, int qual_bin_size):
    """Un-compress the read.
    """
    if not Read_IsCompressed(read):
        return

    if read.seq != NULL:
        uncompress_seq(read, refseq + (read.pos - refstart))

    if read.qual != NULL:
        uncompress_qual(read)

    Read_SetUnCompressed(read)
    return
