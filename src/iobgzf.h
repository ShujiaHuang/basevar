/**
 * @file bgzf.h
 * 
 * @brief C++ wrapper for the codes of htslib/bgzf.h
 * Functions that read and write block gzipped files.
 * 
 * @author Shujia Huang
 * @date 2021-08-20
 * 
 */

#ifndef __INCLUDE_NGSLIB_IOBGZF_H__
#define __INCLUDE_NGSLIB_IOBGZF_H__

#include <iostream>
#include <string>

#include <htslib/bgzf.h>


namespace ngslib {

    class BGZFile {
    private:
        std::string _fname;  // input file name
        std::string _mode;   // read and write mode, Mode matching.
        BGZF *_bgzf;

        /**
         * @brief call `bgzf_open` function to open file.
         * 
         * @param fn 
         * 
         * @param mode  Mode
         * The mode argument can be any of 'r', 'rb', 'a', 'ab', 'w', 'wb', 'x', or 'xb' depending
         * on whether the file will be read or written.  
         * 
         * The default is the mode of fileobj if discernible; otherwise, the default is 'rb'.
         * A mode of 'r' is equivalent to one of 'rb', and similarly for 'w' and 'wb', 'a' and 
         * 'ab', and 'x' and 'xb'.
         * 
         */
        void _open(const std::string &fn, const std::string mode);

        BGZFile(const BGZFile &b) = delete;             // reject using copy constructor (C++11 style).
        BGZFile &operator=(const BGZFile &b) = delete;  // reject using copy/assignment operator (C++11 style).

    public:
        BGZFile(): _bgzf(NULL) {}  // default constructor, do nothing
        explicit BGZFile(const std::string &fn, const std::string mode = "rb") : _fname(fn), _mode(mode) {
            /* Open the specified file for reading or writing. */
            _open(fn, mode);
        }
        ~BGZFile();
        
    };  // class BGZFile
}  // namespace ngslib

#endif

