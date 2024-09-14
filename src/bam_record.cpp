#include <stdexcept>
#include <sstream>

#include <htslib/hts.h>
#include "bam_record.h"
#include "utils.h"


namespace ngslib {

    // _p_cigar_field member should be initialization to a NULL pointer in constructor function.
    BamRecord::BamRecord(const BamRecord &b) : _p_cigar_field(NULL), _n_cigar_op(0) {
        this->_b = bam_dup1(b._b);
        this->_make_cigar_field();
    }

    BamRecord::BamRecord(const bam1_t *b) : _p_cigar_field(NULL), _n_cigar_op(0) {
        this->_b = bam_dup1(b);
        this->_make_cigar_field();
    }

    BamRecord &BamRecord::operator=(const BamRecord &b) {

        if (this->_b)
            bam_destroy1(this->_b);

        this->_b = bam_dup1(b._b);
        this->_make_cigar_field();

        return *this;
    }

    BamRecord &BamRecord::operator=(const bam1_t *b) {

        if (this->_b)
            bam_destroy1(this->_b);

        this->_b = bam_dup1(b);
        this->_make_cigar_field();

        return *this;
    }

    void BamRecord::_make_cigar_field() {
        // 这个函数和 get_cigar_blocks 有相似之处，但我没想好要不要整合
        if (_p_cigar_field)
            delete [] _p_cigar_field;

        if (!_b)
            return;

        _n_cigar_op = _b->core.n_cigar;
        _p_cigar_field = new CigarField[_n_cigar_op];

        if (!_p_cigar_field) {
            throw std::runtime_error("[BamRecord::_make_cigar_field]: Fail to "
                                     "alloc memory space for CigarField.");
        }

        uint32_t *c = bam_get_cigar(_b);
        for (size_t i = 0; i < _n_cigar_op; i++) {
            _p_cigar_field[i].op = bam_cigar_opchr(c[i]);
            _p_cigar_field[i].len = bam_cigar_oplen(c[i]);
        }

        return;
    }

    std::string BamRecord::cigar() const {

        if (!is_mapped()) return "*";  // empty

        std::stringstream cig;
        for (size_t i = 0; i < _n_cigar_op; ++i) {
            cig << _p_cigar_field[i].len << _p_cigar_field[i].op;
        }

        return cig.str();
    }

    void BamRecord::init() {
        if (_b) destroy();
        _b = bam_init1();

        if (_p_cigar_field) {
            delete [] _p_cigar_field;
        }

        _n_cigar_op = 0;
        _p_cigar_field = NULL;

        return;
    }

    int BamRecord::load_read(samFile *fp, sam_hdr_t *h) {

        if (!this->_b)
            this->init();

        int io_status = sam_read1(fp, h, this->_b);

        // Destroy BamRecord and set br to be NULL if fail to read data or
        // hit the end of file.
        if (io_status < 0)
            this->destroy();

        this->_make_cigar_field();
        return io_status;
    }

    int BamRecord::next_read(samFile *fp, hts_itr_t *itr) {

        if (!this->_b)
            this->init();

        int io_status = sam_itr_next(fp, itr, this->_b);

        // Destroy BamRecord and set _b to be NULL if fail to read data or
        // hit the end of file.
        if (io_status < 0)
            this->destroy();

        this->_make_cigar_field();
        return io_status;
    }

    unsigned int BamRecord::align_length() const {

        if (!is_mapped()) return 0;

        unsigned int length = 0;
        char op;
        for (size_t i = 0; i < _n_cigar_op; ++i) {
            op = _p_cigar_field[i].op;
            if (op == 'M' || op == '=' || op == 'X') {
                length += _p_cigar_field[i].len;
            }
        }

        return length;
    }

    unsigned int BamRecord::match_length() const {

        if (!is_mapped()) return 0;

        unsigned int m_size = 0;
        for (size_t i = 0; i < _n_cigar_op; i++) {
            if (_p_cigar_field[i].op == 'M')
                m_size += _p_cigar_field[i].len;
        }

        return m_size;
    }

    void BamRecord::destroy() {
        bam_destroy1(_b);
        _b = NULL;

        if (_p_cigar_field) {
            delete [] _p_cigar_field;
            _p_cigar_field = NULL;
            _n_cigar_op = 0;
        }

        return;
    }

    std::vector<std::tuple<int, uint32_t>> BamRecord::get_cigar_blocks() const {

        if (!_b) {
            throw std::runtime_error("[BamRecord::get_cigar_blocks]: Not found alignment data.");
        }

        std::vector<std::tuple<int, uint32_t>> cigar_block;

        uint32_t *c = bam_get_cigar(_b);
        for (size_t i(0); i < _b->core.n_cigar; ++i) {
            cigar_block.push_back(std::make_tuple(bam_cigar_op(c[i]), bam_cigar_oplen(c[i]))); // (op, l) in sam.h
        }

        return cigar_block;
    }

    std::vector<std::tuple<hts_pos_t, hts_pos_t>> BamRecord::get_alignment_blocks() const {

        // get alignment block by CIGAR
        if (!_b) {
            throw std::runtime_error("[BamRecord::get_alignment_blocks]: Not found alignement data.");
        }

        std::vector<std::tuple<hts_pos_t, hts_pos_t>> align_block;
        hts_pos_t rpos = map_ref_start_pos();  // mapping reference position
        
        int op;         /* cigar OP */ 
        hts_pos_t len;  /* cigar 'OP' length */

        uint32_t *c = bam_get_cigar(_b);
        for (size_t i(0); i < _b->core.n_cigar; ++i) {

            op   = bam_cigar_op(c[i]);
            len  = bam_cigar_oplen(c[i]);

            // 'BAM_XXX' is macros, which defined in sam.h
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {  
                // start position is 0-based
                align_block.push_back(std::make_tuple(rpos, rpos + len));  
                rpos += len;
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                rpos += len;
            }
        }

        return align_block;  // start position is 0-based
    }

    std::vector<ReadAlignedPair> BamRecord::get_aligned_pairs(const std::string &fa) const {
        if (!_b) {
            throw std::runtime_error("[BamRecord::get_aligned_pairs]: Not found alignement data.");
        }

        // mapped pair of read mapped to reference information: 
        // ReadAlignedPair: (cigar_op, read_pos, ref_pos, read_base, read_qual, ref_base)
        std::vector<ReadAlignedPair> aligned_pairs;
        ReadAlignedPair al_pair;

        hts_pos_t rpos = map_ref_start_pos();  // mapping reference position (string index), 0-based
        uint32_t  qpos = 0;                    // mapping read's position (string index), 0-based
        std::string read_seq  = query_sequence();  // read bases
        std::string read_qual = query_qual();      // read bases' qualities

        int op;         /* cigar OP */ 
        hts_pos_t len;  /* cigar 'OP' length */

        uint32_t *c = bam_get_cigar(_b);
        for (size_t i(0); i < _b->core.n_cigar; ++i) {

            op   = bam_cigar_op(c[i]);
            len  = bam_cigar_oplen(c[i]);

            // 'BAM_XXX' is macros, which defined in sam.h
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {  

                for (hts_pos_t i(rpos); i < rpos + len; ++i) {
                    al_pair.op        = op;                         // cigar op
                    al_pair.ref_pos   = i;                          // reference position, 0-based
                    al_pair.ref_base  = fa.substr(i, 1);            // reference base
                    al_pair.qpos      = qpos;                       // read position, 0-based
                    al_pair.read_base = read_seq.substr(qpos, 1);   // read base
                    al_pair.read_qual = read_qual.substr(qpos, 1);  // read quality base

                    aligned_pairs.push_back(al_pair);
                    ++qpos;
                }
                rpos += len;
            } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CPAD) {
                al_pair.op        = op;                           // cigar op
                al_pair.ref_pos   = rpos;                         // reference position, 0-based
                al_pair.ref_base  = "";                           // reference base, empty
                al_pair.qpos      = qpos;                         // read position, 0-based
                al_pair.read_base = read_seq.substr(qpos, len);   // read base
                al_pair.read_qual = read_qual.substr(qpos, len);  // read quality base
                
                aligned_pairs.push_back(al_pair);
                qpos += len;
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                al_pair.op        = op;                    // cigar op
                al_pair.ref_pos   = rpos;                  // reference position, 0-based
                al_pair.ref_base  = fa.substr(rpos, len);  // reference base
                al_pair.qpos      = qpos;                  // read position, 0-based
                al_pair.read_base = "";                    // read base, empty
                al_pair.read_qual = "";   

                aligned_pairs.push_back(al_pair);
                rpos += len;
            } else if (op == BAM_CHARD_CLIP) { 
                /* Nothing to do and get nothing. skipping */ 
                continue;
            }
        }

        return aligned_pairs;  // A vector of tuple.
    }

    unsigned int BamRecord::_max_cigar_Opsize(const char op) const {

        if (!is_mapped()) return 0;

        int max_size = 0;
        for (size_t i = 0; i < _n_cigar_op; i++) {
            if (_p_cigar_field[i].op == op)
                max_size = std::max(_p_cigar_field[i].len, max_size);
        }

        return max_size;
    }

    unsigned int BamRecord::max_insertion_size() const {
        return _max_cigar_Opsize('I');
    }

    unsigned int BamRecord::max_deletion_size() const {
        return _max_cigar_Opsize('D');
    }

    std::string BamRecord::query_sequence() const {

        if (!_b) return "";

        uint8_t *p = bam_get_seq(_b);
        std::string seq(_b->core.l_qseq, 'N');  // initial by a batch of 'N'
        for (size_t i = 0; i < _b->core.l_qseq; ++i)
            seq[i] = _BASES[bam_seqi(p, i)];

        return seq;
    }

    std::string BamRecord::query_qual(int offset) const {

        if (!_b) return "";

        uint8_t *p = bam_get_qual(_b);
        if (!p) return "";

        std::string qual(_b->core.l_qseq, (char)(offset));
        for (size_t i = 0; i < _b->core.l_qseq; ++i)
            qual[i] = (char) (p[i] + offset);

        return qual;
    }

    double BamRecord::mean_qqual() const {

        if (!is_mapped() || (_b->core.l_qseq <= 0))
            return -1;

        double total_phred_score = 0;
        uint8_t *p = bam_get_qual(_b);
        for (size_t i = 0; i < _b->core.l_qseq; ++i)
            total_phred_score += p[i];

        return total_phred_score / _b->core.l_qseq;
    }

    int32_t BamRecord::query_start_pos() const {

        if (!is_mapped()) return -1;

        int32_t p = 0;
        // loop to the end
        for (size_t i = 0; i < _n_cigar_op; ++i) {
            if (_p_cigar_field[i].op == 'S') {
                p += _p_cigar_field[i].len;
            } else {
                break;
            }
        }
        return p;
    }

    int32_t BamRecord::query_end_pos() const {

        if (!is_mapped()) return -1;

        int32_t p = 0;
        // loop from the end
        for (int32_t i = _n_cigar_op - 1; i >= 0; --i) {
            if (_p_cigar_field[i].op == 'S') {
                p += _p_cigar_field[i].len;
            } else { // not a clip, so stop counting
                break;
            }
        }

        return (_b->core.l_qseq - p);
    }

    int32_t BamRecord::query_start_pos_reverse() const {

        if (!is_mapped()) return -1;

        int32_t p = 0;
        // loop from the end
        for (size_t i = _n_cigar_op - 1; i >= 0; --i) {
            if (_p_cigar_field[i].op == 'S') {
                p += _p_cigar_field[i].len;
            } else { // not a clip, stop counting
                break;
            }
        }
        return p;
    }
    
    int32_t BamRecord::query_end_pos_reverse() const {

        if (!is_mapped()) return -1;

        int32_t p = 0;
        for (size_t i = 0; i < _n_cigar_op; ++i) {
            if (_p_cigar_field[i].op == 'S') {
                p += _p_cigar_field[i].len;
            } else { // not a clip, so stop counting
                break;
            }
        }
        return (_b->core.l_qseq - p);
    }

    bool BamRecord::is_proper_orientation() const {

        // _b is NULL
        if (!is_mapped() || !is_mate_mapped()) return false;

        // Get false if mate read mapped on different chromosome
        if (_b->core.tid != _b->core.mtid) return false;

        // Get true if FR: Read1 must map in front of read2, and map
        // to different strands.
        if (_b->core.pos < _b->core.mpos) {
            // Present read is mapped in front of meta => read1.
            // Return true if read1 is mapped to the forward strand (+) and
            // the mate (read2) is mapped to the reverse one (-).
            return ((_b->core.flag & BAM_FREVERSE) == 0) &&
                   ((_b->core.flag & BAM_FMREVERSE) != 0) ? true : false;
        } else {
            // Present read is mapped behind meta => read2
            // Return false if read2 is mapped to forward strand and the
            // mate (read1) is mapped to the reverse one.
            return ((_b->core.flag & BAM_FREVERSE) == 0) &&
                   ((_b->core.flag & BAM_FMREVERSE) != 0) ? false : true;
        }
    }

    bool BamRecord::has_tag(const std::string tag) const {

        if (!_b) return false;

        uint8_t *p = bam_aux_get(_b, tag.c_str());
        return bool(p);
    }

    std::string BamRecord::get_Z_tag(const std::string tag) const {

        std::string tag_str;
        if (has_tag(tag)) {

            uint8_t *p = bam_aux_get(_b, tag.c_str());
            if (*p == 'Z') {
                char *pp = bam_aux2Z(p);
                if (pp) {
                    tag_str = std::string(pp);
                }
            }
        }

        return tag_str;
    }

    std::string BamRecord::get_Int_tag(const std::string tag) const {

        std::string tag_str;
        if (has_tag(tag)) {

            uint8_t *p = bam_aux_get(_b, tag.c_str());
            if (*p == 'I' || *p == 'i' ||
                *p == 'S' || *p == 's' ||
                *p == 'C' || *p == 'c') {

                tag_str = tostring(bam_aux2i(p));
            }
        }

        return tag_str;
    }

    std::string BamRecord::get_Float_tag(const std::string tag) const {

        std::string tag_str;
        if (has_tag(tag)) {

            uint8_t *p = bam_aux_get(_b, tag.c_str());
            if (*p == 'f' || *p == 'd') {
                tag_str = tostring(bam_aux2f(p));
            }
        }

        return tag_str;
    }

    std::string BamRecord::get_tag(const std::string tag) const {

        std::string tag_str;

        tag_str = get_Z_tag(tag);
        if (tag_str.size()) {
            return tag_str;
        }

        tag_str = get_Int_tag(tag);
        if (tag_str.size()) {
            return tag_str;
        }

        tag_str = get_Float_tag(tag);
        if (tag_str.size()) {
            return tag_str;
        }

        return tag_str;
    }

    std::string BamRecord::read_group() const {

        std::string rg;
        if (has_tag("RG")) {
            // try to get from RG tag first
            rg = get_Z_tag("RG");

        } else {
            // try to get the read group tag from qname.
            std::string qn = qname();
            size_t pos = qn.find(":", 0);
            rg = (pos != std::string::npos) ? qn.substr(0, pos) : "";
        }

        // Get the read group, return empty string if no read group found.
        return rg;
    }

    std::ostream &operator<<(std::ostream &os, const BamRecord &r) {

        if (!r._b)
            return os;

        os << r.qname() << "\t"
           << r.flag()  << "\t"
           << r.tid()   << "\t"

           // mapping position +1 to make 1-base coordinate.
           << r.map_ref_start_pos() + 1 << "\t"
           << r.mapq() << "\t"
           << r.cigar() << "\t"
           << r.mate_tid() << "\t"

           // mapping position +1 to make 1-base coordinate.
           << r.mate_map_ref_start_pos() + 1 << "\t"
           << r.insert_size() << "\t"
           << r.query_sequence() << "\t"
           << r.query_qual();

        return os;
    }


}  // namespace ngslib
