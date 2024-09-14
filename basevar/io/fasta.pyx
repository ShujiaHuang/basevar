# cython: profile=True
"""
FastaFile is a utility class used for reading the Fasta file format,
and facilitating access to reference sequences.
"""
import sys
from libc.stdlib cimport atol

from basevar.log import logger
from basevar.io.openfile import Open


cdef class SequenceTuple:
    """Structure for storing data from line of fasta index file.
    """
    cdef public bytes seq_name
    cdef public long int seq_length
    cdef public long int start_position
    cdef public long int line_length
    cdef public long int full_line_length

    def __init__(self, bytes seq_name, long int seq_length, long int start_position,
                 long int line_length, long int full_line_length):
        """Constructor
        """
        self.seq_name = seq_name
        self.seq_length = seq_length
        self.start_position = start_position
        self.line_length = line_length
        self.full_line_length = full_line_length

cdef class FastaIndex:
    """
    Index file of a FastaFile. Contains start and end positions of all
    sequences in the fasta file.

    # .fai format
    # contig, size, location, basesPerLine, bytesPerLine
    # e.g.
    # 1       249250621       52      60      61
    # 2       243199373       253404903       60      61
    """
    def __init__(self, filename, mode="rb", is_ncbi=True):
        """Constructor.

        Takes the name of an Index file. For now, the whole file
        will be read into memory and the references stored in a set.
        """
        self.references = {}
        self.target_name = {}
        self.target_length = {}
        self.n_targets = 0

        self._load_index(filename, mode, is_ncbi)

    def _load_index(self, filename, mode, is_ncbi):

        cdef bytes seq_name
        with Open(filename, mode) as F:

            for the_line in F:

                col = the_line.strip().split()
                seq_name = col[0].split()[0]

                if seq_name.startswith('gi|') and is_ncbi:  # NCBI-formatted line
                    ids = seq_name.split('|')
                    if len(ids) >= 4 and ids[2] == "ref":
                        seq_name = ids[3]

                self.references[seq_name] = SequenceTuple(col[0], atol(col[1]), atol(col[2]), atol(col[3]), atol(col[4]))

                self.target_name[self.n_targets] = self.references[seq_name].seq_name
                self.target_length[self.n_targets] = self.references[seq_name].seq_length

                # n_targets is the number of reference sequences
                self.n_targets += 1
        return


cdef class FastaFile:
    """
    Utility for reading sequence from Fasta files.
    """
    def __init__(self, fastafile, indexfile, mode="rb", parseNCBI=True):
        """
        Constructor. Takes file-name and index file-name
        """
        self.filename = fastafile
        self.the_file = Open(fastafile, mode)
        self.the_index = FastaIndex(indexfile, mode=mode, is_ncbi=parseNCBI)

        self.references = self.the_index.references
        self.cache_ref_name = None
        self.cache_start_pos = -1
        self.cache_end_pos = -1
        self.cache = None

    def get_total_sequence_length(self):
        """
        Return the accumulated lengths of all sequences in the
        reference index.
        """
        return sum([seq_tuple.seq_length for _, seq_tuple in self.references.iteritems()])

    def get_reference_length(self, seq_name):

        try:
            return self.references[seq_name].seq_length
        except KeyError:
            raise KeyError, "Cannot find %s in reference fasta file" % seq_name

    cpdef void close(self):
        """
        Wrapper function to close self.theFile
        """
        self.the_file.close()

    cdef char *get_character(self, bytes seq_name, long int pos):
        """
        Returns the character at the specified (0-indexed means 0-base system) position
        of the specified sequence.
        """
        if self.cache is not None:
            if self.cache_start_pos <= pos < self.cache_end_pos:
                # logger.debug("Getting %s:%s-%s from cache. cache index = %s:%s" % (
                #     seq_name, pos, pos, pos - self.cache_start_pos, pos - self.cache_start_pos))

                return self.cache[pos - self.cache_start_pos]

        cdef SequenceTuple seq_tuple = self.references[seq_name]
        cdef long int seq_length = seq_tuple.seq_length
        cdef long int seq_start_position = seq_tuple.start_position
        cdef long int line_length = seq_tuple.line_length
        cdef long int full_line_length = seq_tuple.full_line_length

        if pos >= seq_length or pos < 0:
            # it's empty
            return <char*> ""

        filepos = seq_start_position + pos + (full_line_length - line_length) * (
            <long int>((<double>pos)/line_length))
        self.the_file.seek(filepos)

        try:
            return self.the_file.read(1)
        except Exception:
            # return nothing
            return <char*> ""

    cdef void set_cache_sequence(self, bytes seq_name, long int begin_pos, long int end_pos):
        """cache a sequence in memery make the program much faster
        """
        if seq_name not in self.references:
            logger.error("Invalid contig name %s. Make sure your FASTA reference file and query regions "
                         "have the same naming convention" % seq_name)
            sys.exit(1)

        cdef SequenceTuple seq_tuple = self.references[seq_name]
        cdef long int seq_length = seq_tuple.seq_length

        # it's 0-base system
        begin_pos = max(0, begin_pos)
        end_pos = min(seq_length - 1, end_pos)

        cdef long int seq_start_pos = seq_tuple.start_position
        cdef long int line_length = seq_tuple.line_length
        cdef long int full_line_length = seq_tuple.full_line_length
        cdef long int desired_sequence_start_pos = seq_start_pos + begin_pos + (
                full_line_length - line_length) * <long int> ((<double> begin_pos) / line_length)
        cdef long int desired_sequence_end_pos = seq_start_pos + end_pos + (
                full_line_length - line_length) * <long int> ((<double> end_pos) / line_length)

        # 0-base
        cdef long int desired_seq_length = (end_pos - begin_pos)
        cdef long int desired_sequence_length_infile = (desired_sequence_end_pos - desired_sequence_start_pos)

        if end_pos < begin_pos:
            raise IndexError, "Cannot have beginPos = %s, endPos = %s" % (begin_pos, end_pos)

        if end_pos > seq_length or begin_pos < 0:
            raise IndexError, ("Cannot return sequence from %s to %s. Reference sequence "
                               "length = %s" % (begin_pos, end_pos, seq_length))

        self.the_file.seek(desired_sequence_start_pos)
        self.cache = self.the_file.read(desired_sequence_length_infile).replace("\n", "")
        self.cache_ref_name = seq_name
        self.cache_start_pos = begin_pos
        self.cache_end_pos = end_pos

    cdef bytes get_sequence(self, bytes seq_name, long int begin_pos, long int end_pos):
        """
        Returns the character sequence between the the specified (0-indexed) start
        and end positions. This returns a half-open sequence interval, i.e. the returned
        sequence includes the character at beginPos, but not the one at endPos. This is done
        in order to make down-stream sequence handling easier.
        """
        if self.cache is not None:
            if begin_pos >= self.cache_start_pos and end_pos < self.cache_end_pos:
                # logger.debug("Getting %s:%s-%s from cache. cache index = %s:%s" % (
                #     seq_name, begin_pos, end_pos, begin_pos - self.cache_start_pos, end_pos - self.cache_start_pos))
                return self.cache[begin_pos - self.cache_start_pos:end_pos - self.cache_start_pos]

        cdef SequenceTuple seq_tuple = self.references[seq_name]
        cdef long int seq_length = seq_tuple.seq_length

        # make 0-base system
        begin_pos = max(0, begin_pos-1)
        end_pos = min(seq_length - 1, end_pos)

        cdef long int seq_start_pos = seq_tuple.start_position
        cdef long int line_length = seq_tuple.line_length
        cdef long int full_line_length = seq_tuple.full_line_length
        cdef long int desired_sequence_start_pos = seq_start_pos + begin_pos + (
                full_line_length - line_length) * <long int> ((<double> begin_pos) / line_length)
        cdef long int desired_sequence_end_pos = seq_start_pos + end_pos + (
                full_line_length - line_length) * <long int> ((<double> end_pos) / line_length)

        # 0-base
        cdef long int desired_seq_length = (end_pos - begin_pos)
        cdef long int desired_sequence_length_infile = (desired_sequence_end_pos - desired_sequence_start_pos)

        if end_pos < begin_pos:
            raise IndexError, "Cannot have beginPos = %s, endPos = %s" % (begin_pos, end_pos)

        if end_pos > seq_length or begin_pos < 0:
            raise IndexError, ("Cannot return sequence from %s to %s. Reference sequence "
                               "length = %s" % (begin_pos, end_pos, seq_length))

        self.the_file.seek(desired_sequence_start_pos)
        cdef bytes seq = self.the_file.read(desired_sequence_length_infile)
        return seq.replace("\n", "")

    property filename:
        """The filename of the `reference` sequences file. """
        def __get__(self):
            return self.filename

    property nreferences:
        """number of `reference` sequences in the file."""
        def __get__(self):
            return self.the_index.n_targets

    property refnames:
        """tuple with the name of the `reference` sequences"""
        def __get__(self):
            s = []
            cdef int x = 0
            for x in range(self.the_index.n_targets):
                s.append(self.the_index.target_name[x])

            return tuple(s)

    property lengths:
        """tuple of the lengths of the `reference` sequences. The length are in the same order as `FastaFile.refnames` 
        """
        def __get__(self):
            s = []
            cdef int x = 0
            for x in range(self.the_index.n_targets):
                s.append(self.the_index.target_length[x])

            return tuple(s)

