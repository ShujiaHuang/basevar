#include <stdexcept>

#include <htslib/hts.h>
#include "bam_header.h"
#include "utils.h"

namespace ngslib {

    BamHeader::BamHeader(samFile *fp) {
        _h = sam_hdr_read(fp);
    }

    BamHeader::BamHeader(const std::string &fn) {

        if (!is_readable(fn)) {
            throw std::runtime_error("_bam_header::BamHeader: " + fn + " not found.");
        }

        samFile *fp = hts_open(fn.c_str(), "r");
        _h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);
    }

    BamHeader &BamHeader::operator=(const BamHeader &bh) {
        // release _h pointer if _h is not NULL
        sam_hdr_destroy(_h);
        _h = sam_hdr_dup(bh._h);
        return *this;
    }

    BamHeader &BamHeader::operator=(const sam_hdr_t *hdr) {
        // release _h pointer if _h is not NULL
        sam_hdr_destroy(_h);
        _h = sam_hdr_dup(hdr);
        return *this;
    }

    std::ostream &operator<<(std::ostream &os, const BamHeader &hd) {

        if (hd._h)
            os << sam_hdr_str(hd._h);

        return os;
    }

    void BamHeader::destroy() {
        sam_hdr_destroy(_h);
        _h = NULL;
    }

    int BamHeader::name2id(const std::string &name) const {
        int tid = sam_hdr_name2tid(_h, name.c_str());

        if (tid < 0) {
            throw std::runtime_error(
                "[bam_header.cpp::BamHeader:name2id] Unknown reference name or "
                "the header not be parsed: " + name);
        }
        return tid;
    }

    std::string BamHeader::get_sample_name() {
        size_t r, n_rg = sam_hdr_count_lines(_h, "RG");
        kstring_t sm = KS_INITIALIZE;  // {0, 0, NULL}
        for (size_t i(0); i < n_rg; ++i) {
            // 0 on success; 
            // -1 if the requested tag does not exist;
            // -2 on other errors
            r = sam_hdr_find_tag_pos(_h, "RG", i, "SM", &sm);
            if (r < 0) continue;  // not found, continue
            break;                // 一个 bam 里可能有多个 RG，但找到第一个 SM 后就退出，不管其他的
        }

        std::string samplename = sm.s ? std::string(sm.s) : "";
        if (sm.s == NULL) {
            throw std::runtime_error("[bam_header.cpp::BamHeader:get_sample_name] Bam file format error: "
                                     "missing `SM` tag in `@RG` field in BAM/CRAM/SAM header.");
        } else {
            free(sm.s);
        }

        return samplename;
    }

}  // namespace ngslib