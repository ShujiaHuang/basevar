# cython: profile=True
"""Wrapper for htslib
"""
import os
from warnings import warn
from libc.errno cimport errno
from posix.unistd cimport dup

from basevar.io.libcutils cimport encode_filename
from basevar.io.libcutils cimport force_str, charptr_to_str


__all__ = ['HTSFile', 'Samfile', 'ReadIterator', 'destroy_read']


# defines imported from samtools
DEF SEEK_SET = 0
DEF SEEK_CUR = 1
DEF SEEK_END = 2

# maximum genomic coordinace
cdef int MAX_POS = (1 << 31) - 1

cdef tuple FORMAT_CATEGORIES = ('UNKNOWN', 'ALIGNMENTS', 'VARIANTS', 'INDEX', 'REGIONS')
cdef tuple FORMATS = ('UNKNOWN', 'BINARY_FORMAT', 'TEXT_FORMAT', 'SAM', 'BAM', 'BAI', 'CRAM', 'CRAI',
                      'VCF', 'BCF', 'CSI', 'GZI', 'TBI', 'BED')
cdef tuple COMPRESSION = ('NONE', 'GZIP', 'BGZF', 'CUSTOM')


# valid types for sam headers
VALID_HEADER_TYPES = {"HD": dict,
                      "SQ": list,
                      "RG": list,
                      "PG": list,
                      "CO": list}


# type conversions within sam header records
cdef dict VALID_HEADER_FIELDS = {
    "HD": {"VN": str, "SO": str, "GO": str},
    "SQ": {"SN": str, "LN": int, "AH": str, "AS": str, "M5": str, "UR": str, "SP": str},
    "RG": {"ID": str, "SM": str, "LB": str, "DS": str, "PU": str, "PI": str, "CN": str, "DT": str,
           "PL": str, "PG": str},
    "PG": {"ID": str, "VN": str, "CL": str, "PN": str, "PP": str, "DS": str}, }


cdef int COMPRESS_COUNT = 40


########################################################################
########################################################################
## Helpers for backward compatibility to hide the difference between
## boolean properties and methods
########################################################################

class CallableValue(object):
    def __init__(self, value):
        self.value = value
    def __call__(self):
        return self.value
    def __bool__(self):
        return self.value
    def __nonzero__(self):
        return self.value
    def __eq__(self, other):
        return self.value == other
    def __ne__(self, other):
        return self.value != other


CTrue = CallableValue(True)
CFalse = CallableValue(False)


cdef class HTSFile(object):
    """
    Base class for HTS file types
    """

    def __cinit__(self, *args, **kwargs):
        self.htsfile = NULL
        self.threads = 1
        self.duplicate_filehandle = True

    def close(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL

    def check_truncation(self, ignore_truncation=False):
        """Check if file is truncated."""
        if not self.htsfile:
            return

        if self.htsfile.format.compression != bgzf:
            return

        cdef BGZF *bgzfp = hts_get_bgzfp(self.htsfile)
        if not bgzfp:
            return

        cdef int ret = bgzf_check_EOF(bgzfp)
        if ret < 0:
            raise IOError(errno, 'error checking for EOF marker')
        elif ret == 0:
            msg = 'no BGZF EOF marker; file may be truncated'.format(self.filename)
            if ignore_truncation:
                warn(msg)
            else:
                raise IOError(msg)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    @property
    def category(self):
        """General file format category.  One of UNKNOWN, ALIGNMENTS,
        VARIANTS, INDEX, REGIONS"""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return FORMAT_CATEGORIES[self.htsfile.format.category]

    @property
    def format(self):
        """File format.

        One of UNKNOWN, BINARY_FORMAT, TEXT_FORMAT, SAM, BAM,
        BAI, CRAM, CRAI, VCF, BCF, CSI, GZI, TBI, BED.
        """
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return FORMATS[self.htsfile.format.format]

    @property
    def version(self):
        """Tuple of file format version numbers (major, minor)"""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return self.htsfile.format.version.major, self.htsfile.format.version.minor

    @property
    def compression(self):
        """File compression.

        One of NONE, GZIP, BGZF, CUSTOM."""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return COMPRESSION[self.htsfile.format.compression]

    @property
    def description(self):
        """Vaguely human readable description of the file format"""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        cdef char *desc = hts_format_description(&self.htsfile.format)
        try:
            return charptr_to_str(desc)
        finally:
            free(desc)

    @property
    def is_open(self):
        """return True if HTSFile is open and in a valid state."""
        return CTrue if self.htsfile != NULL else CFalse

    @property
    def is_closed(self):
        """return True if HTSFile is closed."""
        return self.htsfile == NULL

    @property
    def closed(self):
        """return True if HTSFile is closed."""
        return self.htsfile == NULL

    @property
    def is_write(self):
        """return True if HTSFile is open for writing"""
        return self.htsfile != NULL and self.htsfile.is_write != 0

    @property
    def is_read(self):
        """return True if HTSFile is open for reading"""
        return self.htsfile != NULL and self.htsfile.is_write == 0

    @property
    def is_sam(self):
        """return True if HTSFile is reading or writing a SAM alignment file"""
        return self.htsfile != NULL and self.htsfile.format.format == sam

    @property
    def is_bam(self):
        """return True if HTSFile is reading or writing a BAM alignment file"""
        return self.htsfile != NULL and self.htsfile.format.format == bam

    @property
    def is_cram(self):
        """return True if HTSFile is reading or writing a BAM alignment file"""
        return self.htsfile != NULL and self.htsfile.format.format == cram

    @property
    def is_vcf(self):
        """return True if HTSFile is reading or writing a VCF variant file"""
        return self.htsfile != NULL and self.htsfile.format.format == vcf

    @property
    def is_bcf(self):
        """return True if HTSFile is reading or writing a BCF variant file"""
        return self.htsfile != NULL and self.htsfile.format.format == bcf

    def reset(self):
        """reset file position to beginning of file just after the header.

        Returns
        -------

        The file position after moving the file pointer.

        """
        return self.seek(self.start_offset)

    def seek(self, uint64_t offset):
        """move file pointer to position *offset*, see :meth:`pysam.HTSFile.tell`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise IOError('seek not available in streams')

        cdef int64_t ret
        if self.htsfile.format.compression == bgzf:
            with nogil:
                ret = bgzf_seek(hts_get_bgzfp(self.htsfile), offset, SEEK_SET)
        elif self.htsfile.format.compression == no_compression:
            with nogil:
                ret = hts_useek(self.htsfile, <int>offset, SEEK_SET)
        else:
            raise NotImplementedError("seek not implemented in files compressed by method {}".format(
                self.htsfile.format.compression))
        return ret

    def tell(self):
        """return current file position, see :meth:`pysam.HTSFile.seek`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise IOError('tell not available in streams')

        cdef int64_t ret
        if self.htsfile.format.compression == bgzf:
            with nogil:
                ret = bgzf_tell(hts_get_bgzfp(self.htsfile))
        elif self.htsfile.format.compression == no_compression:
            with nogil:
                ret = hts_utell(self.htsfile)
        elif self.htsfile.format.format == cram:
            with nogil:
                ret = htell(cram_fd_get_fp(self.htsfile.fp.cram))
        else:
            raise NotImplementedError("seek not implemented in files compressed by method {}".format(
                self.htsfile.format.compression))

        return ret

    cdef htsFile *_open_htsfile(self) except? NULL:
        cdef char *cfilename
        cdef char *cmode = self.mode
        cdef int fd, dup_fd, threads

        threads = self.threads - 1
        if isinstance(self.filename, bytes):
            cfilename = self.filename
            with nogil:
                htsfile = hts_open(cfilename, cmode)
                if htsfile != NULL:
                    hts_set_threads(htsfile, threads)
                return htsfile
        else:
            if isinstance(self.filename, int):
                fd = self.filename
            else:
                fd = self.filename.fileno()

            if self.duplicate_filehandle:
                dup_fd = dup(fd)
            else:
                dup_fd = fd

            # Replicate mode normalization done in hts_open_format
            smode = self.mode.replace(b'b', b'').replace(b'c', b'')
            if b'b' in self.mode:
                smode += b'b'
            elif b'c' in self.mode:
                smode += b'c'
            cmode = smode

            hfile = hdopen(dup_fd, cmode)
            if hfile == NULL:
                raise IOError('Cannot create hfile')

            try:
                # filename.name can be an int
                filename = str(self.filename.name)
            except AttributeError:
                filename = '<fd:{}>'.format(fd)

            filename = encode_filename(filename)
            cfilename = filename
            with nogil:
                htsfile = hts_hopen(hfile, cfilename, cmode)
                if htsfile != NULL:
                    hts_set_threads(htsfile, threads)
                return htsfile

    def add_hts_options(self, format_options=None):
        """Given a list of key=value format option strings, add them to an open htsFile
        """
        cdef int rval
        cdef hts_opt *opts = NULL

        if format_options:
            for format_option in format_options:
                rval = hts_opt_add(&opts, format_option)
                if rval != 0:
                    if opts != NULL:
                        hts_opt_free(opts)
                    raise RuntimeError('Invalid format option ({}) specified'.format(format_option))
            if opts != NULL:
                rval = hts_opt_apply(self.htsfile, opts)
                if rval != 0:
                    hts_opt_free(opts)
                    raise RuntimeError('An error occured while applying the requested format options')
                hts_opt_free(opts)

    def parse_region(self, contig=None, start=None, stop=None,
                     region=None, tid=None,
                     reference=None, end=None):
        """parse alternative ways to specify a genomic region. A region can
        either be specified by :term:`contig`, `start` and
        `stop`. `start` and `stop` denote 0-based, half-open
        intervals. :term:`reference` and `end` are also accepted for
        backward compatiblity as synonyms for :term:`contig` and
        `stop`, respectively.

        Alternatively, a samtools :term:`region` string can be
        supplied.

        If any of the coordinates are missing they will be replaced by
        the minimum (`start`) or maximum (`stop`) coordinate.

        Note that region strings are 1-based inclusive, while `start`
        and `stop` denote an interval in 0-based, half-open
        coordinates (like BED files and Python slices).

        If `contig` or `region` or are ``*``, unmapped reads at the end
        of a BAM file will be returned. Setting either to ``.`` will
        iterate from the beginning of the file.

        Returns
        -------

        tuple : a tuple of `flag`, :term:`tid`, `start` and
        `stop`. The flag indicates whether no coordinates were
        supplied and the genomic region is the complete genomic space.

        Raises
        ------

        ValueError
           for invalid or out of bounds regions.

        """
        cdef int rtid
        cdef int32_t rstart
        cdef int32_t rstop

        if reference is not None:
            if contig is not None:
                raise ValueError('contig and reference should not both be specified')
            contig = reference

        if end is not None:
            if stop is not None:
                raise ValueError('stop and end should not both be specified')
            stop = end

        if contig is None and tid is None and region is None:
            return 0, 0, 0, MAX_POS

        rtid = -1
        rstart = 0
        rstop = MAX_POS
        if start is not None:
            try:
                rstart = start
            except OverflowError:
                raise ValueError('start out of range (%i)' % start)

        if stop is not None:
            try:
                rstop = stop
            except OverflowError:
                raise ValueError('stop out of range (%i)' % stop)

        if region:
            region = force_str(region)
            if ":" in region:
                contig, coord = region.split(":")
                parts = coord.split("-")
                rstart = int(parts[0]) - 1
                if len(parts) >= 1:
                    rstop = int(parts[1])
            else:
                contig = region

        if tid is not None:
            if not self.is_valid_tid(tid):
                raise IndexError('invalid tid')
            rtid = tid
        else:
            if contig == "*":
                rtid = HTS_IDX_NOCOOR
            elif contig == ".":
                rtid = HTS_IDX_START
            else:
                rtid = self.get_tid(contig)
                if rtid < 0:
                    raise ValueError('invalid contig `%s`' % contig)

        if rstart > rstop:
            raise ValueError('invalid coordinates: start (%i) > stop (%i)' % (rstart, rstop))
        if not 0 <= rstart < MAX_POS:
            raise ValueError('start out of range (%i)' % rstart)
        if not 0 <= rstop <= MAX_POS:
            raise ValueError('stop out of range (%i)' % rstop)

        return 1, rtid, rstart, rstop

    def is_valid_tid(self, tid):
        """
        return True if the numerical :term:`tid` is valid; False otherwise.

        returns -1 if contig is not known.
        """
        raise NotImplementedError()

    def is_valid_reference_name(self, contig):
        """
        return True if the contig name :term:`contig` is valid; False otherwise.
        """
        return self.get_tid(contig) != -1

    def get_tid(self, contig):
        """
        return the numerical :term:`tid` corresponding to
        :term:`contig`

        returns -1 if contig is not known.
        """
        raise NotImplementedError()

    def get_reference_name(self, tid):
        """
        return :term:`contig` name corresponding to numerical :term:`tid`
        """
        raise NotImplementedError()


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
    def __cinit__(self, char* filename):
        """Constructor.
        """
        self.filename = filename
        self.samfile = NULL
        self.the_header = NULL
        self.index = NULL

    def __dealloc__(self):
        """Clean up.
        Here I've pasted the code from several clean-up methods, as calling methods
        is not guaranteed to work in this function.
        """
        # remember: dealloc cannot call other methods
        if self.samfile != NULL:
            sam_close(self.samfile)
            self.samfile = NULL

        self.clear_index()
        self.clear_header()

    cdef void clear_header(self):
        """ Clear all the header data. """
        if self.the_header != NULL:
            bam_hdr_destroy(self.the_header)
            self.the_header = NULL

    cdef void clear_index(self):
        """Clear all the index file data. """
        if self.index != NULL:
            hts_idx_destroy(self.index)
            self.index = NULL

    cdef bint _is_bam(self):
        return self.samfile.is_bin

    cdef bint _is_cram(self):
        return self.samfile.is_cram

    cdef bint _is_open(self):
        """return true if samfile has been opened."""
        return self.samfile != NULL

    cdef bint _has_index(self):
        """return true if samfile has an existing (and opened) index."""
        return self.index != NULL

    cdef void open(self, basestring mode, bint load_index):
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
                                  "Check that index exists." % self.filename)

    cdef void _open_bamfile(self, mode):
        """Open BamFile.
        """
        self.samfile = sam_open(self.filename, mode)
        self.the_header = bam_hdr_init()
        self.the_header = sam_hdr_read(self.samfile)

        return

    cdef char* getrname(self, int tid):
        """Convert numerical :term:`tid` into :ref:`reference` name."""
        if not 0 <= tid < self.the_header.n_targets:
            raise ValueError("tid (%s) out of range 0<=tid<%i" % (tid, self.the_header.n_targets))

        return self.the_header.target_name[tid]

    cdef ReadIterator fetch(self, const char *region):
        """Fetch reads from a specified region."""

        # Need to load the index for new queries.
        if not self._is_open():
            self.open("r", True)

        if self._is_bam() or self._is_cram():
            return ReadIterator(self, region)
        else:
            raise StandardError, "Random access query only allowed for BAM/CRAM files."

    cdef void close(self):
        """closes file."""
        if not self._is_cram():
            if self.samfile != NULL:
                sam_close(self.samfile)
                self.samfile = NULL

        self.clear_index()
        self.clear_header()

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
        """tuple of the lengths of the :term:`reference` sequences. 
        The lengths are in the same order as :attr:`pysam.Samfile.reference`.
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

        ``region`` is looks like: chr:start-end
        """

        self.the_samfile = NULL
        self.the_iterator = NULL
        self.b = NULL

        if not samfile._is_open():
            raise StandardError, "Samfile %s is not open. Cannot read from file." % (samfile.filename)

        self.the_samfile = samfile.samfile

        # Load from BAM by querying index for file off-set
        if samfile._has_index():
            self.the_iterator = sam_itr_querys(samfile.index, samfile.the_header, region)
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
        """ Some very importance data structure in sam.h and wrapper by htslibWrapper module.
        
        /*************************
         *** Alignment records ***
         *************************/

        /*! @typedef
         @abstract Structure for core alignment information.
         @field  tid     chromosome ID, defined by bam_hdr_t
         @field  pos     0-based leftmost coordinate
         @field  bin     bin calculated by bam_reg2bin()
         @field  qual    mapping quality
         @field  l_qname length of the query name
         @field  flag    bitwise flag
         @field  l_extranul length of extra NULs between qname & cigar (for alignment)
         @field  n_cigar number of CIGAR operations
         @field  l_qseq  length of the query sequence (read)
         @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
         @field  mpos    0-based leftmost coordinate of next read in template
         */
        typedef struct {
            int32_t tid;
            int32_t pos;
            uint16_t bin;
            uint8_t qual;
            uint8_t l_qname;
            uint16_t flag;
            uint8_t unused1;
            uint8_t l_extranul;
            uint32_t n_cigar;
            int32_t l_qseq;
            int32_t mtid;
            int32_t mpos;
            int32_t isize;
        } bam1_core_t;
        
        Notes: read pos is 0-base!
        """
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
            seq[i] = self._get_base(s, i)
            qual[i] = q[i]
            assert qual[i] <= 93
            assert qual[i] >= 0

        seq[len_seq] = '\0'
        qual[len_seq] = '\0'

        read_start = c.pos  # 0-base
        cdef short* cigar_ops = <short*> malloc(2 * c.n_cigar * sizeof(short))
        assert cigar_ops != NULL, "Error cigar_ops is NULL"

        cdef uint32_t* cigar = bam_get_cigar(self.b)
        for i in range(c.n_cigar):
            cigar_flag = bam_cigar_op(cigar[i])
            cigar_flag_len = bam_cigar_oplen(cigar[i])

            cigar_ops[2 * i] = cigar_flag
            cigar_ops[(2 * i) + 1] = cigar_flag_len

            # Soft-clipping of sequence at start of read changes the mapping
            # position. Recorded mapping pos is that of the first aligned (not soft-clipped)
            # base. I want to adjust this so that the read start refers to the first base in
            # the read.
            if i == 0 and cigar_flag == 4:
                read_start -= cigar_flag_len

        the_read.seq = seq
        the_read.qual = qual
        the_read.cigar_ops = cigar_ops
        the_read.hash = NULL
        the_read.mate_chrom_id = c.mtid
        the_read.cigar_len = c.n_cigar
        the_read.chrom_id = c.tid
        the_read.r_len = len_seq
        the_read.pos = read_start
        the_read.end = bam_endpos(self.b)
        the_read.insert_size = c.isize
        the_read.mate_pos = c.mpos
        the_read.bit_flag = c.flag
        the_read.mapq = c.qual

        Read_SetUnCompressed(the_read)

        return the_read

    cdef int cnext(self) nogil:
        """cversion of iterator. Used by IteratorColumn."""
        return sam_itr_next(self.the_samfile, self.the_iterator, self.b) >= 0

    cdef char _get_base(self, uint8_t *s, int i):
        cdef char* base_lookup = "=ACMGRSVTWYHKDBN"
        return base_lookup[bam_seqi(s, i)]


cdef void destroy_read(cAlignedRead* the_read):
    """De-allocate memory for read.
    """

    if the_read.seq != NULL:
        free(the_read.seq)

    if the_read.qual != NULL:
        free(the_read.qual)

    if the_read.cigar_ops != NULL:
        free(the_read.cigar_ops)

    if the_read.hash != NULL:
        free(the_read.hash)

    free(the_read)


cdef void compress_seq(cAlignedRead* read, char* refseq):
    """Does exactly what it says on the tin.
    """
    cdef char* seq = read.seq
    cdef char* new_seq = <char*>(calloc(read.r_len * 2, sizeof(char)))
    cdef char* final_seq = NULL
    cdef int i = 0
    cdef int n_matches = 0
    cdef int new_seq_index = 0
    cdef int total_matches = 0

    for i in range(read.r_len):

        # Ref match. Store count if > COMPRESS_COUNT
        if seq[i] == refseq[i]:
            total_matches += 1

            if n_matches == COMPRESS_COUNT:
                new_seq[new_seq_index] = n_matches
                n_matches = 0
                new_seq_index += 1

            n_matches += 1

        # Ref mis-match. Store base.
        else:
            if n_matches > 0:
                new_seq[new_seq_index] = n_matches
                n_matches = 0
                new_seq_index += 1

            new_seq[new_seq_index] = seq[i]
            new_seq_index += 1

    # If we finish off with a load of matches, catch those
    if n_matches > 0:
        new_seq[new_seq_index] = n_matches
        n_matches = 0
        new_seq_index += 1

    new_seq[new_seq_index] = 0
    final_seq = <char*>calloc(new_seq_index + 1, sizeof(char))
    strcpy(final_seq, new_seq)

    final_seq[new_seq_index] = 0
    free(read.seq)
    free(new_seq)

    read.seq = final_seq
    return


cdef void compress_qual(cAlignedRead* read, int qual_bin_size):
    """Does exactly what it says on the tin.
    """
    cdef char* qual = read.qual
    cdef char* new_qual = <char*>(calloc(2 * read.r_len, sizeof(char)))
    cdef char* final_qual = NULL
    cdef char last_char = 0
    cdef int last_count = 0
    cdef int new_qual_index = 0

    cdef int i = 0
    if qual_bin_size > 1:
        for i in range(read.r_len):
            qual[i] = (qual[i] / qual_bin_size) * qual_bin_size

    for i in range(read.r_len):

        # assert qual[i] >= 0, "Shit 3. Qual = %s" % (qual[i])
        # assert qual[i] <= 93, "Shit 4!. Qual = %s" % (qual[i])

        if i == 0:
            new_qual[new_qual_index] = qual[i] + 33
            new_qual_index += 1
            last_char = qual[i]
            last_count = 1
        else:
            if qual[i] == last_char:
                last_count += 1
            else:
                new_qual[new_qual_index] = last_count
                new_qual_index += 1
                new_qual[new_qual_index] = qual[i] + 33
                new_qual_index += 1
                last_char = qual[i]
                last_count = 1

    if last_count > 0:
        new_qual[new_qual_index] = last_count
        new_qual_index += 1

    final_qual = <char*>malloc(sizeof(char) * (new_qual_index + 1))
    new_qual[new_qual_index] = 0
    strcpy(final_qual, new_qual)

    final_qual[new_qual_index] = 0
    free(read.qual)
    free(new_qual)

    read.qual = final_qual
    return


cdef void uncompress_seq(cAlignedRead* read, char* refseq):
    """Does exactly what it says on the tin.
    """
    cdef int ref_index = 0
    cdef int new_seq_index = 0
    cdef int i = 0
    cdef int j = 0
    cdef int seq_len = strlen(read.seq)
    cdef char* seq = read.seq
    cdef char* new_seq = <char*>(malloc(sizeof(char*) * (read.r_len + 1)))

    for i in range(seq_len):
        if seq[i] <= COMPRESS_COUNT:
            for j in range(seq[i]):
                new_seq[new_seq_index] = refseq[j + ref_index]
                new_seq_index += 1

            ref_index += seq[i]

        else:
            new_seq[new_seq_index] = seq[i]
            ref_index += 1
            new_seq_index += 1

    free(read.seq)
    read.seq = new_seq
    read.seq[read.r_len] = 0
    return


cdef void uncompress_qual(cAlignedRead* read):
    """
    Does exactly what it says on the tin.
    """
    cdef char* new_qual = <char*>(malloc(sizeof(char) * (read.r_len + 1)))
    cdef char* qual = read.qual

    cdef int len_qual = strlen(read.qual)
    cdef int i = 0
    cdef int j = 0
    cdef int new_qual_index = 0

    for i in range(0, len_qual - 1, 2):
        for j in range(qual[i + 1]):
            new_qual[new_qual_index] = qual[i] - 33
            assert qual[i] - 33 >= 0, "Shit 1. Qual = %s. Count = %s!" % (qual[i], qual[i + 1])
            assert qual[i] - 33 <= 93, "Shit 2!. Qual = %s. Count = %s" % (qual[i], qual[i + 1])
            new_qual_index += 1

    assert read.r_len == new_qual_index

    free(read.qual)
    read.qual = new_qual
    read.qual[read.r_len] = 0
    return


cdef void compress_read(cAlignedRead* read, char* refseq, int refstart, int refend, int qual_bin_size):
    """To save memory use, we compress reads. The sequence is compressed using reference-based compression.
    
    The qualities are compressed using zlib and (optionally) using course binning. We delete the read-group tag, 
    and the cigar string information (this can only be done after the read has been used for candidate generation).
    A bit is set in the bit-field to label the read as compressed or uncompressed.
    """
    if Read_IsCompressed(read):
        return

    if read.seq != NULL:
        compress_seq(read, refseq + (read.pos - refstart))

    if read.qual != NULL:
        compress_qual(read, qual_bin_size)

    if read.hash != NULL:
        free(read.hash)
        read.hash = NULL

    Read_SetCompressed(read)
    return


cdef void uncompress_read(cAlignedRead* read, char* refseq, int refstart, int refend, int qual_bin_size):
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
