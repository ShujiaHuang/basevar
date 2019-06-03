"""A Wrapper for htslib

Author: Shujia Huang
Date: 2019-05-29 18:20:35
"""
from cpython cimport bool

cdef extern from "stdlib.h":
    void *malloc(size_t)
    void *calloc(size_t, size_t)
    void free(void *)
    void *alloca(size_t)
    int c_abs "abs"(int)

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)

cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass

    FILE *fopen(char *, char *)
    FILE *freopen(char *path, char *mode, FILE *stream)
    int fileno(FILE *stream)
    int dup2(int oldfd, int newfd)
    int fflush(FILE *stream)
    FILE *stderr
    FILE *stdout
    int fclose(FILE *)
    int sscanf(char *str, char *fmt, ...)
    int printf(char *str, char *fmt, ...)
    int sprintf(char *str, char *fmt, ...)
    int fprintf(FILE *ifile, char *fmt, ...)
    char *fgets(char *str, int size, FILE *ifile)

cdef extern from "ctype.h":
    int toupper(int c)

cdef extern from "string.h":
    ctypedef int size_t
    void *memcpy(void *dst, void *src, size_t len)
    int strncmp(char *s1, char *s2, size_t len)
    char *strcpy(char *dest, char *src)
    char *strncpy(char *dest, char *src, size_t len)
    size_t strlen(char *s)
    int memcmp(void *s1, void *s2, size_t len)

cdef extern from "stdint.h":
    ctypedef int int8_t
    ctypedef int int64_t
    ctypedef int int32_t
    ctypedef int uint32_t
    ctypedef int uint8_t
    ctypedef int uint64_t

cdef extern from "htslib/bgzf.h":
    ctypedef struct z_stream:
        pass

    ctypedef struct hFILE:
        pass

    ctypedef struct bgzf_mtaux_t:
        pass

    ctypedef struct bgzidx_t:
        pass

    ctypedef struct BGZF:
        int errcode
        int is_write
        int is_be
        int compress_level
        int is_compressed
        int is_gzip
        int cache_size
        int block_length
        int block_offset
        int64_t block_address
        int64_t uncompressed_address
        void *uncompressed_block
        void *compressed_block
        void *cache
        hFILE *fp
        bgzf_mtaux_t *mt
        bgzidx_t *idx
        int idx_build_otf
        z_stream *gz_stream

cdef extern from "htslib/hts.h":
    ctypedef struct kstring_t:
        size_t l
        size_t m
        char *s

    ctypedef struct hts_idx_t:
        pass

    ctypedef struct hts_itr_t:
        pass

    void hts_idx_destroy(hts_idx_t *idx)
    const char *hts_parse_reg(const char *s, int *beg, int *end)

cdef extern from "htslib/sam.h":
    ctypedef struct bam_hdr_t:
        int32_t n_targets
        int32_t ignore_sam_err
        uint32_t l_text
        uint32_t *target_len
        int8_t *cigar_tab
        char **target_name
        char *text
        void *sdict

    ctypedef struct bam1_core_t:
        int32_t tid
        int32_t pos
        uint32_t bin
        uint32_t qual
        uint32_t l_qname
        uint32_t flag
        uint32_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int32_t mpos
        int32_t isize

    ctypedef struct bam1_t:
        bam1_core_t core
        int l_data
        int m_data
        uint8_t *data
        uint64_t id

    ctypedef struct cram_fd:
        pass

    ctypedef union samFileUnion:
        BGZF *bgzf
        cram_fd *cram
        hFILE *hfile
        void *voidp

    ctypedef struct samFile:
        uint32_t is_bin
        uint32_t is_write
        uint32_t is_be
        uint32_t is_cram
        uint32_t is_compressed
        uint32_t is_kstream
        uint32_t dummy
        int64_t lineno
        kstring_t line
        char *fn
        char *fn_aux
        samFileUnion fp

    samFile *sam_open(const char *fn, const char *mode)
    int bam_name2id(bam_hdr_t *h, const char *ref)
    hts_idx_t *sam_index_load(samFile *fp, const char *fn)
    bam_hdr_t *bam_hdr_init()
    bam_hdr_t *sam_hdr_read(samFile *fp)
    bam1_t *bam_init1()
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b)
    bool bam_is_rev(const bam1_t *b)
    bool bam_is_mrev(const bam1_t *b)
    char *bam_get_qname(const bam1_t *b)
    uint8_t *bam_get_seq(const bam1_t *b)
    uint8_t *bam_get_qual(const bam1_t *b)
    uint8_t *bam_get_aux(const bam1_t *b)
    uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
    char *bam_aux2Z(const uint8_t *s)
    uint32_t *bam_get_cigar(const bam1_t *b)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
    int sam_itr_next(samFile *htsfp, hts_itr_t *itr, bam1_t *r) nogil
    int sam_close(samFile *fp)
    uint8_t bam_seqi(uint8_t *s, int i)
    void bam_destroy1(bam1_t *b)
    void bam_hdr_destroy(bam_hdr_t *h)
    void sam_itr_destroy(hts_itr_t *itr)
    int32_t bam_endpos(const bam1_t *b)
    int bam_cigar_op(uint32_t c)
    int bam_cigar_oplen(uint32_t c)


ctypedef struct cAlignedRead:
    char *seq
    char *qual
    short *cigar_ops
    short *hash
    short mate_chrom_id
    short cigar_len
    short chrom_id
    short r_len
    int pos
    int end
    int insert_size
    int mate_pos
    int bit_flag
    unsigned char mapq


cdef class ReadIterator:
    cdef cAlignedRead *get(self, int store_rgID, char** rgID)
    cdef int cnext(self) nogil
    cdef char _get_base(self, uint8_t *s, int i)

    cdef samFile *the_samfile
    cdef hts_itr_t *the_iterator
    cdef bam1_t *b


cdef class Samfile:
    cdef void clear_header(self)
    cdef void clear_index(self)

    cpdef void _open(self, mode, bool load_index)  # A main function to open BAM/CRAM/SAM
    cdef void _open_bamfile(self, mode)
    cdef int _is_bam(self)
    cdef int _is_cram(self)
    cdef int _is_open(self)
    cdef int _has_index(self)  # return 1 (for True) or 0 (for False)

    cdef char* getrname(self, int tid)
    cdef ReadIterator fetch(self, const char *region)
    cpdef close(self)

    cdef char* filename
    cdef samFile *samfile
    cdef bam_hdr_t *the_header
    cdef hts_idx_t *index
    cdef object lock


###################################################################################################
# These are bits set in the BAM bit-flag.
DEF BAM_FPAIRED = 1  # Read is paired in sequencing, no matter whether it is mapped in a pair
DEF BAM_FPROPER_PAIR = 2  # Read is mapped in a proper pair
DEF BAM_FUNMAP = 4  # Read itself is unmapped; conflictive with BAM_FPROPER_PAIR
DEF BAM_FMUNMAP = 8  # Mate is unmapped
DEF BAM_FREVERSE = 16  # Read is mapped to the reverse strand
DEF BAM_FMREVERSE = 32  # Mate is mapped to the reverse strand
DEF BAM_FREAD1 = 64  # This is read1
DEF BAM_FREAD2 = 128  # This is read2
DEF BAM_FSECONDARY = 256  # Not primary alignment
DEF BAM_FQCFAIL = 512  # QC failure
DEF BAM_FDUP = 1024  # Optical or PCR duplicate
DEF BAM_FCOMPRESSED = 2048  # Is the read compressed
###################################################################################################


# And here are accessor functions for the bit-fields
cdef inline int Read_IsReverse(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FREVERSE) != 0)

cdef inline int Read_IsPaired(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FPAIRED) != 0)

cdef inline int Read_IsProperPair(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FPROPER_PAIR) != 0)

cdef inline int Read_IsDuplicate(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FDUP) != 0)

cdef inline int Read_IsUnmapped(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FUNMAP) != 0)

cdef inline int Read_MateIsUnmapped(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FMUNMAP) != 0)

cdef inline int Read_MateIsReverse(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FMREVERSE) != 0)

cdef inline int Read_IsQCFail(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FQCFAIL) != 0)

cdef inline int Read_IsReadOne(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FREAD1) != 0)

cdef inline int Read_IsSecondaryAlignment(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FSECONDARY) != 0)

cdef inline int Read_IsCompressed(cAlignedRead* the_read) nogil:
    return ((the_read.bit_flag & BAM_FCOMPRESSED) != 0)

cdef inline int Read_SetIsNotReverse(cAlignedRead* the_read) nogil:
    the_read.bit_flag &= (~BAM_FREVERSE)

cdef inline int Read_SetIsReverse(cAlignedRead* the_read) nogil:
    the_read.bit_flag |= BAM_FREVERSE

cdef inline void Read_SetQCFail(cAlignedRead* the_read) nogil:
    the_read.bit_flag |= BAM_FQCFAIL

cdef inline void Read_SetCompressed(cAlignedRead* the_read) nogil:
    the_read.bit_flag |= BAM_FCOMPRESSED

cdef inline void Read_SetUnCompressed(cAlignedRead* the_read) nogil:
    the_read.bit_flag &= (~BAM_FCOMPRESSED)

###################################################################################################

cdef void destroy_read(cAlignedRead* the_read)
cdef void compress_read(cAlignedRead* read, char* refSeq, int ref_start, int refend, int qual_bin_size, int full_comp)
cdef void uncompress_read(cAlignedRead* read, char* refSeq, int ref_start, int refend, int qual_bin_size)
