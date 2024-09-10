// The C++ codes for dealing BAM header
// Author: Shujia Huang
// Date: 2021-08-20

#ifndef __INCLUDE_NGSLIB_BAM_HEADER_H__
#define __INCLUDE_NGSLIB_BAM_HEADER_H__

#include <iostream>
#include <string>

#include <htslib/sam.h>


namespace ngslib {

    // Store the header of BAM file, which also acts as a dictionary of
    // reference sequences with names and lengths.
    class BamHeader {

    private:

        /*! `sam_hdr_t` is defined in sam.h. The data structure of `sam_hdr_t` is:
         *
         * @abstract Structure for the alignment header.
         *
         * @field n_targets   number of reference sequences
         * @field l_text      length of the plain text in the header (may be zero if
         *                    the header has been edited)
         * @field target_len  lengths of the reference sequences
         * @field target_name names of the reference sequences
         * @field text        plain text (may be NULL if the header has been edited)
         * @field sdict       header dictionary
         * @field hrecs       pointer to the extended header struct (internal use only)
         * @field ref_count   reference count
         * @note The text and l_text fields are included for backwards compatibility.

         These fields may be set to NULL and zero respectively as a side-effect
         of calling some header API functions.  New code that needs to access the
         header text should use the sam_hdr_str() and sam_hdr_length() functions
         instead of these fields.
         */
        sam_hdr_t *_h;   // `bam_hdr_t` is an old name of `sam_hdr_t`, do not use it.

    public:

        // Initializes a new empty BamHeader with no data.
        BamHeader() : _h(NULL) {}
        ~BamHeader() { destroy(); }

        // Only explicit conversions allowed.
        explicit BamHeader(samFile *fp);

        /** Read the header from a BAM compressed file.
         * This function works on SAM, BAM and CRAM files.
         */
        explicit BamHeader(const std::string &fn);

        // Create BamHeader from a exist header, rarely use.
        BamHeader(const sam_hdr_t *hdr) { _h = sam_hdr_dup(hdr); }

        // Copy constructor.
        BamHeader(const BamHeader &bh) { _h = sam_hdr_dup(bh._h); }
        BamHeader &operator=(const BamHeader &bh);
        BamHeader &operator=(const sam_hdr_t *hdr);
        friend std::ostream &operator<<(std::ostream &os, const BamHeader &hd);

        void init() {
            if (_h) destroy();
            _h = sam_hdr_init();
        }

        // Free the memory and set Bam file header pointer to be NULL (could save memory).
        void destroy();

        operator bool() const { return bool(_h != NULL); }

        // Write BAM header to a BAM file.
        int write(samFile *fp) {
            // samFile is an alias of htsFile which define in sam.h
            return sam_hdr_write(fp, _h);
        }

        // return the `sam_hdr_t` pointer of BAM file header.
        sam_hdr_t *h() const { return _h; }

        // Return the names of the reference sequences by the index of chromosome
        // in header.
        std::string seq_name(int i) const {
            return std::string(_h->target_name[i]);
        }

        /// Get the target id for a given reference sequence name
        /*!
         * @param ref  Reference name
         * @return     Positive value on success,
         *             -1 if unknown reference,
         *             -2 if the header could not be parsed
         *
         * Looks up a reference sequence by name in the reference hash table
         * and returns the numerical target id.
         */
        int name2id(const std::string &name) const;

        // Return a length of the reference sequences by the index of chromosome
        uint32_t seq_length(int i) const { return _h->target_len[i]; }

        // Returns the text representation of the header.
        std::string header_txt() {
            return std::string(sam_hdr_str(_h));
        }

        // Get sample ID from SM tag in @RG field.
        std::string get_sample_name();
    };
}

#endif