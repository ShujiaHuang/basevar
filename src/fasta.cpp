#include <stdexcept>

#include "fasta.h"
#include "utils.h"


namespace ngslib {

    void Fasta::_load_data(const char *file_name) {
        fai = NULL;

        if (!is_readable(file_name)) {
            throw std::invalid_argument("fasta::Fasta: file not found - " + tostring(file_name));
        }

        fname = tostring(file_name);
        fai = fai_load(file_name);  // load data

        if (!fai) {
            throw std::invalid_argument("fasta::Fasta: index not loaded.");
        }
    }

    Fasta::Fasta(const Fasta &ft) {  // copy constructor
        this->_load_data(ft.fname.c_str());   // re-load FASTA file.
    }

    Fasta &Fasta::operator=(const char *file_name) {
        if (fai) {
            fai_destroy(fai);
        }

        this->_load_data(file_name);  // re-load the file
        return *this;
    }

    // return the sequence string of seq_id
    std::string &Fasta::operator[](const std::string reg) {

        std::map<std::string, std::string>::iterator it;
        it = _seq.find(reg);
        if (it == _seq.end()) {
            // record the seq into _seq if not 
            _seq[reg] = fetch(reg);
        }

        return _seq[reg];
    }

    std::string Fasta::fetch(const char *reg) const {
        // check if we have loaded the fasta index
        if (!fai) throw std::invalid_argument("Fasta::fetch index not loaded");

        int length;
        char *f = fai_fetch(fai, reg, &length);

        if (!f) {
            throw std::invalid_argument("Fasta::fetch - Fail to fetch sequence.");
        }

        std::string sub_seq(f);
        free(f);
        
        if (sub_seq.empty()) {
            throw std::invalid_argument("Fasta::fetch - Fetch empty sequence on " + tostring(reg));
        }

        return sub_seq;
    }

    std::string Fasta::fetch(const char *chromosome,
                             const uint32_t start,   // start: 0-base, end: 0-base.
                             const uint32_t end) const {

        // check if we have loaded the fasta index
        if (!fai) throw std::invalid_argument("Fasta::fetch index not loaded");
        if (start > end) throw std::invalid_argument("Fasta::fetch the start position must be <= end.");
        if (start < 0) throw std::invalid_argument("Fasta::fetch the start position must be >= 0");

        int length;
        char *f = faidx_fetch_seq(fai, chromosome, start, end, &length);

        if (!f) {
            throw std::invalid_argument("Fasta::fetch - Fail to fetch sequence.");
        }

        std::string sub_seq(f);
        free(f);

        if (sub_seq.empty()) {
            throw std::invalid_argument("Fasta::fetch - Fetch empty sequence on " + tostring(chromosome) +
                                        ":" + tostring(start) + "-" + tostring(end));
        }
        return sub_seq;
    }

    // Output the filename of FASTA
    std::ostream &Fasta::len_out(std::ostream &os) const {
        if (fai) {
            os << "\nThe length information: \n";

            int nseq = faidx_nseq(fai);
            const char *seq_name;
            for (int i = 0; i < nseq; i++) {
                seq_name = faidx_iseq(fai, i);
                os << seq_name << " = " << faidx_seq_len(fai, seq_name) << "\n";
            }
        }

        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Fasta &fa) {
        os << fa.fname;
        fa.len_out(os);  // use private method for output
        return os;
    }

}  // namespace ngslib


