#include "iobgzf.h"
#include "utils.h"

namespace ngslib {
    void BGZFile::_open(const std::string &fn, const std::string mode) {

        if (fn.find('t') != std::string::npos || fn.find('U') != std::string::npos)
            throw std::invalid_argument("[iobgzf.cpp::BGZFile:_open] Invalid mode: " + mode);
        if ((mode.find('r') && !is_readable(fn)))
            throw std::runtime_error("[iobgzf.cpp::BGZFile:_open] file not found - " + fn);
        
        _bgzf = bgzf_open(fn.c_str(), mode.c_str());
        if (!_bgzf) {
            throw std::runtime_error("[iobgzf.cpp::BGZFile:_open] file open failure.");
        }

        return;
    }

    BGZFile::~BGZFile() {
        if (_bgzf) bgzf_close(_bgzf);  // `bgzf_close` return 0 on success and -1 on error
    }
    
}  // namespace ngslib