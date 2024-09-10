// C++ codes for dealing BAM data
// Author: Shujia Huang
// Date: 2021-08-22

#ifndef __INCLUDE_NGSLIB_BAM_RECORD_H__
#define __INCLUDE_NGSLIB_BAM_RECORD_H__

#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include <htslib/sam.h>
#include "bam_header.h"


namespace ngslib {

    /*! _BASES is defined according to information of `bam_get_seq()` in sam.h,
     *  detail for bam_get_seq() is bellow:
     *
     *  @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
     *  8 for T and 15 for N. Two bases are packed in one byte with the base
     *  at the higher 4 bits having smaller coordinate on the read. It is
     *  recommended to use bam_seqi() macro to get the base.
     *
     **/
    static const char _BASES[16] = {' ', 'A', 'C', ' ',
                                    'G', ' ', ' ', ' ',
                                    'T', ' ', ' ', ' ',
                                    ' ', ' ', ' ', 'N'};

    /// Defined a structure for recording the mapped pair information of read mapped on reference.
    // been used in `BamRecord::get_aligned_pairs` 
    typedef struct {
        int         op;         // cigar op,    detail about: `bam_get_cigar(_b)` could be found in sam.h
        hts_pos_t   ref_pos;    // reference position
        std::string ref_base;   // reference base
        uint32_t    qpos;       // read position
        std::string read_base;  // read base
        std::string read_qual;  // read quality base
    } ReadAlignedPair;             

    class BamRecord {

    private:
        /// bam record, the most and only important member of BamRecord
        bam1_t *_b;

        /** Define a structure for CIGAR op by type BAM_CIGAR_STR (MIDNSHP=XB)
         * and length.
         * @param op    CIGAR op (MIDNSHP=XB)
         * @param len   CIGAR length
         * */
        typedef struct {
            char op;
            int len;
        } CigarField;

        // A pointer to CigarField array.
        CigarField *_p_cigar_field;

        // The number of CIGAR operator, which is the size of CigarField array.
        unsigned int _n_cigar_op;

        /* Make cigar field by CIGAR of this alignment */
        void _make_cigar_field();

        /* get the max size of Op in cigar */
        unsigned int _max_cigar_Opsize(const char op) const;

    public:
        BamRecord() : _b(NULL), _p_cigar_field(NULL), _n_cigar_op(0) {}  // default constructor, initial to be NULL.
        ~BamRecord() { destroy(); }

        BamRecord(const BamRecord &b);  // copy constructor
        BamRecord &operator=(const BamRecord &b);

        BamRecord(const bam1_t *b);
        BamRecord &operator=(const bam1_t *b);

        void init();
        void destroy();

        /**************************
         *** Exported functions ***
         **************************/

        /// The two TOP-level functions communicate with outside source file.
        /// _b will be changed only in these two functions.

        /** Call sam_read1() function of sam.h to load a record from file.
         *
         *  @param fp   Pointer to the source file
         *  @param h    Pointer to the header previously read (fully or partially)
         *  @return >= 0 on successfully reading a new record, -1 on end of
         *  stream, < -1 on error.
         *
         **/
        int load_read(samFile *fp, sam_hdr_t *h);

        /** Call sam_itr_next() function of sam.h to get iterator the next
         *  read by a SAM/BAM/CRAM iterator.
         *
         *  @param fp       samFile pointer to the source file
         *  @param itr      Iterator
         *  @return >= 0 on success; -1 when there is no more data; < -1 on error.
         *
         **/
        int next_read(samFile *fp, hts_itr_t *itr);

        // -- END the two TOP level functions --

        // conversion function
        operator bool() const { return bool(_b != NULL); }
        friend std::ostream &operator<<(std::ostream &os, const BamRecord &b);

        /// 12 inline functions for dealing with FLAG of BAM alignment record

        /* Reads are Pair-end in sequencing, no matter whether it is
         * mapped in a pair.
         * */
        bool is_paired() const { return _b && (_b->core.flag & BAM_FPAIRED); }

        /* This is read1 */
        bool is_read1() const { return _b && (_b->core.flag & BAM_FREAD1); }

        /* This is read2 */
        bool is_read2() const { return _b && (_b->core.flag & BAM_FREAD2); }

        /* The read itself is mapped */
        bool is_mapped() const { return _b && ((_b->core.flag & BAM_FUNMAP) == 0); }

        /* The read is mapped to the reverse strand (-) */
        bool is_mapped_reverse() const {
            return is_mapped() && (_b->core.flag & BAM_FREVERSE);
        }

        /* The meta read is mapped */
        bool is_mate_mapped() const {
            return is_paired() && ((_b->core.flag & BAM_FMUNMAP) == 0);
        }

        /* The mate read is mapped to the reverse strand (-) */
        bool is_mate_mapped_reverse() const {
            return is_mate_mapped() && (_b->core.flag & BAM_FMREVERSE);
        }

        /*  The read is mapped in a proper pair */
        bool is_proper_pair() const {
            return is_paired() && (_b->core.flag & BAM_FPROPER_PAIR);
        }

        /* The read is a secondary alignment (not primary) */
        bool is_secondary() const {
            return is_mapped() && (_b->core.flag & BAM_FSECONDARY);
        }

        /* The mapped read is failed QC */
        bool is_qc_fail() const {
            return is_mapped() && (_b->core.flag & BAM_FQCFAIL);
        }

        /* The read is a duplicate */
        bool is_duplicate() const { return is_mapped() && (_b->core.flag & BAM_FDUP); }

        /* The read is a supplementary alignment */
        bool is_supplementary() const {
            return is_mapped() && (_b->core.flag & BAM_FSUPPLEMENTARY);
        }

        // -- End the FLAG inline functions --

        /// Functions need input a Bam header

        /* Get the alignment chromosome of this read */
        std::string tid_name(const BamHeader &hdr) const {
            return is_mapped() ? hdr.seq_name(_b->core.tid) : "";
        }

        // hts_pos_t is a alisa name of int64_t defined in hts.h.
        hts_pos_t tid_length(const BamHeader &hdr) const {
            return is_mapped() ? hdr.seq_length(_b->core.tid) : -1;
        }

        /* Get the alignment chromosome of mate read */
        std::string mate_tid_name(const BamHeader &hdr) const {
            return is_mate_mapped() ? hdr.seq_name(_b->core.mtid) : "";
        }

        hts_pos_t mate_tid_length(const BamHeader &hdr) const {
            return is_mate_mapped() ? hdr.seq_length(_b->core.mtid) : -1;
        }

        /// Functions for alignment reference information.

        /* Get the full alignment flag of this read */
        uint16_t flag() const { return _b->core.flag; }

        /* Get the id of alignment chromosome, defined by sam_hdr_t */
        int32_t tid() const { return is_mapped() ? _b->core.tid : -1; }

        /* chromosome ID of mate read in template, defined by sam_hdr_t */
        int32_t mate_tid() const { return is_mate_mapped() ? _b->core.mtid : -1; }

        /* Get the alignment strand of this read, Should be one of '*', '-', '+' */
        char map_strand() const {
            return is_mapped() ? (is_mapped_reverse() ? '-' : '+') : '*';
        }

        /* Get the alignment strand of mate read, Should be one of '*', '-', '+' */
        char mate_map_strand() const {
            return is_mate_mapped() ? (is_mate_mapped_reverse() ? '-' : '+') : '*';
        }

        /* Get the begin mapped position on the reference genome, 0-based coordinate.
         * @return  a hts_pos_t value on success, -1 on NULL.
         * */
        hts_pos_t map_ref_start_pos() const {
            return is_mapped() ? _b->core.pos : -1;
        }

        /* Get the end mapped position of the alignment on the reference.
         * Calculate the rightmost base position of an alignment on the reference
         * genome.
         *
         * @return  The coordinate of the first base after the alignment, 1-based.
         *
         * For a mapped read, this is just b->core.pos + bam_cigar2rlen.
         * For an unmapped read (either according to its flags or if it has no cigar
         * string) or a read whose cigar string consumes no reference bases at all,
         * we return b->core.pos + 1 by convention.
         */
        hts_pos_t map_ref_end_pos() const {
            return is_mapped() ? bam_endpos(_b) : -1;
        }

        /* Get the begin mapped position of mate on the reference genome, 0-based
         * @return  a hts_pos_t value on success, -1 on NULL.
         * */
        hts_pos_t mate_map_ref_start_pos() const {
            // 0-based leftmost coordinate of next read in template
            return is_mate_mapped() ? _b->core.mpos : -1;
        }

        /* Get mapping quality */
        int mapq() const { return is_mapped() ? _b->core.qual : 0; }

        /* convert CIGAR to a string */
        std::string cigar() const;

        /*
         * Return the number of "aligned bases" exclude the base mark as I, D, N,
         * S, H, and P in CIGAR.
         */
        unsigned int align_length() const;

        /* Get the total number of matched base ('M') in this alignment. */
        unsigned int match_length() const;

        /* Get max insertion size of this alignment. */
        unsigned int max_insertion_size() const;

        /* Get max deletion size of this alignment. */
        unsigned int max_deletion_size() const;

        /* Get insert size */
        hts_pos_t insert_size() const { return is_paired() ? _b->core.isize : 0; }

        /// Functions for the alignment query information
        /* Get the read ID as a string */
        std::string qname() const { return std::string(bam_get_qname(_b)); }

        /* Get the length of read */
        int query_length() const { return _b ? _b->core.l_qseq : -1; }

        /* Retrieve the sequencing bases of this read as a string (ACTGN) */
        std::string query_sequence() const;

        /* Retrieve the sequencing qualities of this read as a string
         *
         * @param offset Encoding offset for Phred quality scores. Default 33
         * @return quality scores after converting offset. If first char is empty,
         * returns empty string.
         *
         * */
        std::string query_qual(int offset = 33) const;

        /* Calculate the mean sequencing quality of the whole read */
        double mean_qqual() const;

        /* Get the alignment start position on this read, by removing soft-clips.
         *
         * @return 0-base position on the read on success, -1 on NULL.
         *
         * */
        int32_t query_start_pos() const;

        /* Get the end of the alignment on the read, by removing soft-clips, 1-base
         * @return The last position in 1-base on the read on success, -1 on NULL.
         * */
        int32_t query_end_pos() const;

        /* Get the alignment start position on this read, by removing soft-clips.
         *
         * Do it in the reverse orientation.
         * @return 0-base position on the read on success, -1 on NULL.
         *
         * */
        int32_t query_start_pos_reverse() const;

        /* Get the end of the alignment on the read, by removing soft-clips, 1-base
         * Do it in the reverse orientation.
         * @return The last position in 1-base on the read on success, -1 on NULL.
         * */
        int32_t query_end_pos_reverse() const;

        /**
         * @brief Get the cigar block object
         * 
         * The cigar string order in the array is "MIDNSHP=X" followed by a field for the NM tag. 
         * If the NM tag is not present, this field will always be 0. 
         * 
         * Detail in htslib/sam.h: "CIGAR related macros".
         * 
            +-----+--------------+-----+
            |M    |BAM_CMATCH    |0    |
            +-----+--------------+-----+
            |I    |BAM_CINS      |1    |
            +-----+--------------+-----+
            |D    |BAM_CDEL      |2    |
            +-----+--------------+-----+
            |N    |BAM_CREF_SKIP |3    |
            +-----+--------------+-----+
            |S    |BAM_CSOFT_CLIP|4    |
            +-----+--------------+-----+
            |H    |BAM_CHARD_CLIP|5    |
            +-----+--------------+-----+
            |P    |BAM_CPAD      |6    |
            +-----+--------------+-----+
            |=    |BAM_CEQUAL    |7    |
            +-----+--------------+-----+
            |X    |BAM_CDIFF     |8    |
            +-----+--------------+-----+
            |B    |BAM_CBACK     |9    |
            +-----+--------------+-----+
            |NM   |NM tag        |10   |
            +-----+--------------+-----+
         * 
         * If the cigar string is : '3M1I50M5D46M'
         * then alignment.cigar is: [(0, 3), (1, 1), (0, 50), (2, 5), (0, 46)]
         * 
         * @return std::vector<std::tuple<int, uint32_t>> 
         * 
         */
        std::vector<std::tuple<int, uint32_t>> get_cigar_blocks() const;

        /**
         * @brief Get a vector of tuple of start and end positions of aligned gapless blocks.

            The start and end positions are in genomic coordinates. 
            The start postion is 0-based, end position is 1-based.

            Blocks are not normalized, i.e. two blocks might be directly adjacent. 
            This happens if the two blocks are separated by an insertion in the read.

            If the cigar string is: 3M1I50M5D46M, and the start mapping position is 22862751 in BAM.
            Alignment.blocks looks like: [(22862750, 22862753), (22862753, 22862803), (22862808, 22862854)]

            Ref the `get_block(self)` function in: https://github.com/pysam-developers/pysam/blob/master/pysam/libcalignedsegment.pyx
        * 
        * @return std::vector<std::tuple<hts_pos_t, hts_pos_t>> 
        * 
        */
        std::vector<std::tuple<hts_pos_t, hts_pos_t>> get_alignment_blocks() const;

        /// mapped pair of read mapped to reference information: (cigar_op, read_pos, ref_pos, read_base, read_qual, ref_base)
        /**
         * @brief Get a vector of aligned read (query) and reference positions. 
         * 
         * Each item in the returned vector is a structure of ReadAlignedPair consisting of the 0-based offset 
         * from the start of the read sequence followed by the 0-based reference position.
         * 
         * Ref the `get_aligned_pairs` function in https://github.com/pysam-developers/pysam/blob/master/pysam/libcalignedsegment.pyx
         * 
         * @param fa 
         * 
         * @return 
         * std::vector<ReadAlignedPair> 
         * (cigar_op, read_pos, ref_pos, read_base, read_qual, ref_base)
         * 
         */
        std::vector<ReadAlignedPair> get_aligned_pairs(const std::string &fa) const;

        /// Other useful functions
        /* BamRecord has proper orientation (FR): lower position read is mapped to
         * forward strand(+) the higher one mapped to reverse strand(-).
         * */
        bool is_proper_orientation() const;

        /* has a specific TAG in the alignment or not */
        bool has_tag(const std::string tag) const;

        /** Get a string (Z) tag
         * @param tag Name of the tag. eg "XP"
         * @param s The string to be filled in with the tag information
         * @return the tag value as a string is present, even if empty.
         * */
        std::string get_Z_tag(const std::string tag) const;

        /** Get an int (i) tag
         * @param tag Name of the tag. eg "XP"
         * @param t Value to be filled in with the tag value.
         * @return the tag value as a string is present, even if empty.
         * */
        std::string get_Int_tag(const std::string tag) const;

        /** Get a float (f) tag
         * @param tag Name of the tag. eg "AS"
         * @param t Value to be filled in with the tag value.
         * @return the tag value as a string is present, even if empty.
         * */
        std::string get_Float_tag(const std::string tag) const;

        /** Get a string of either Z, f or i type. Useful if tag type not known
         * at compile time.
         *
         * @param tag Name of the tag. eg "XP"
         * @param s The string to be filled in with the tag information
         * @return the tag value as a string is present, even if empty.
         *
         * */
        std::string get_tag(const std::string tag) const;

        /** Get the read group, first from RG tag , then by qname.
         * @return empty string if no read group found
         */
        std::string read_group() const;

        /// Set data to the alignment to change the status of alignment record
        /* Set QC fail for this alignment read */
        void set_qc_fail() {
            if (is_mapped()) _b->core.flag |= BAM_FQCFAIL;
        }

    };  // class BamRecord

}  // namespace ngslib

#endif
