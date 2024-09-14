// The C++ code for reading FASTA.
// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_FASTA_H__
#define __INCLUDE_NGSLIB_FASTA_H__

#include <iostream>
#include <string>
#include <map>

#include <htslib/faidx.h>

namespace ngslib {

    // Identify the FASTA format
    class Fasta {
    private:
        std::string fname;
        faidx_t *fai;
        std::map<std::string, std::string> _seq;

        // Load the FASTA indexed of reference sequence. The index file (.fai) will be build if
        // the reference file doesn't have one. The input file could be bgzip-compressed.
        void _load_data(const char *name);
        std::ostream &len_out(std::ostream &os) const;

    public:
        // default constructor
        Fasta() : fai(NULL) {}
        Fasta(const char *file_name) { this->_load_data(file_name); }
        Fasta(const std::string &file_name) { this->_load_data(file_name.c_str()); }
        Fasta(const Fasta &);  // copy constructor

        // Destroy the malloc'ed faidx_t index inside object
        ~Fasta() { if (fai) fai_destroy(fai); }

        Fasta &operator=(const char *s);
        Fasta &operator=(const std::string &s) { return *this = s.c_str(); }  // inline definition
        Fasta &operator=(const Fasta &s) { return *this = s.fname; }  // inline definition

        // Return the sequence string of reg in format "chr2:2-10"
        std::string &operator[](const std::string reg);
        friend std::ostream &operator<<(std::ostream &os, const Fasta &fa);

        // Query if sequence is present
        /* @param  fai  Pointer to the faidx_t struct
         * @param  seq  Sequence name
         * @return      1 if present or 0 if absent
         */
        bool has_seq(const std::string seq_id) const { return faidx_has_seq(fai, seq_id.c_str()); }

        // Return number of sequences in fa index
        int nseq() { return faidx_nseq(fai); }

        // Return name of i-th sequence
        std::string iseq_name(int i) { return std::string(faidx_iseq(fai, i)); } 

        // Return sequence length, -1 if not present
        uint32_t seq_length(const char *seq_id) const { 
            return faidx_seq_len64(fai, seq_id) > 0 ? faidx_seq_len64(fai, seq_id) : 0; 
        }
        uint32_t seq_length(const std::string seq_id) const { return seq_length(seq_id.c_str()); }

        /** 
         * @brief fetch fasta sequence in a region.
         * 
         * @param reg Region in the format "chr2:20,000-30,000", "chr2:20,000" or "chr2"  (1-based)
         * @return std::string 
         * 
         * Note: The start position in `reg` must be 1-based, and the real fetch is [start-1, end).
         * 
         */
        std::string fetch(const char *reg) const;
        std::string fetch(const std::string &reg) const { return fetch(reg.c_str()); }

        /** fetch a string from the fasta sequence
         * @param chromosome name of the reference to query
         * @param start position.  1. Zero-based
         * @param end position.    2. Zero-based
         *
         * @exception Throws an invalid_argument if start > end, chromosome not found, or seq not found
         * @note This is currently NOT thread safe
         */
        std::string fetch(const char *chromosome, const uint32_t start, const uint32_t end) const;
        std::string fetch(const std::string &chromosome, const uint32_t start, const uint32_t end) const {
            return fetch(chromosome.c_str(), start, end);
        }
        std::string fetch(const std::string &chromosome, const uint32_t start) const {
            return fetch(chromosome, start, seq_length(chromosome));
        }
    };  // class Fasta

}  // namespace ngslib

#endif  // #ifndef __INCLUDE_NGSLIB_FASTA_H__
