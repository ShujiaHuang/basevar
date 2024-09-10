// Author: Shujia Huang
// Date: 2021-08-22
#include <iostream>
#include <string>

#include "bam_header.h"

int main() {

    using ngslib::BamHeader;
    std::string fn1 = "../data/range.bam";
    std::string fn2 = "../data/range.cram";
    std::string fn3 = "../data/xx_minimal.sam";

    BamHeader bh0;
    BamHeader bh1(fn1);
    // BamHeader bh2 = fn2; // Not allow
    BamHeader bh2 = BamHeader("../data/range.cram");

    // sam_hdr_t assign to BamHeader.
    BamHeader bh3(bh2.h());
    BamHeader bh4;
    bh4 = bh3;
    bh4 = bh2.h();

    // Test NULL to NULL
    BamHeader bh5;
    BamHeader bh6 = bh5;
    bh6 = bh5;

    // bh7 and bh8 is equally.
    BamHeader bh7 = BamHeader(fn3);
    BamHeader bh8 = BamHeader(fn3);

    std::cout << ">>> The header of bh0 is NULL: " << bh0 << "\n";
    if (!bh0) {
        std::cout << "bh0 is empty: " << !bh0 << "; bool: " << bool(bh0) << "\n";
    }

    std::cout << ">>> 1. The header of file: " + fn1 << "\n" << bh1 << "\n";
    std::cout << ">>> 2. The header of file: " + fn2 << "\n" << bh2 << "\n";
    std::cout << ">>> 3. The header of file: bh3\n" << bh3 << "\n";
    std::cout << ">>> 4. The header of file: bh4\n" << bh4 << "\n";
    std::cout << ">>> 5. The header of file: bh5" << bh5 << "; bool: " << bool(bh5) << "\n";
    std::cout << ">>> 6. The header of file: bh6" << bh6 << "; bool: " << bool(bh6) << "\n";
    std::cout << ">>> 7. The header of file: bh7\n" << bh7 << "\n";
    std::cout << ">>> 8. The header of file: bh8\n" << bh8 << "\n";

    // samFile assign to BamHeader
    samFile *fp = sam_open(fn1.c_str(), "r");
    bh8 = BamHeader(fp);  // Now the file pointer has moved to other place.
    std::cout << ">>> 9. The header of file: bh8\n" << bh8 << "\n";
    BamHeader bh9 = bh8;
    bh8.destroy();
    std::cout << ">>> 10. The header of file: bh9\n" << bh9 << "\n";

    std::cout << "bh9.seq_name(1): " << bh9.seq_name(1)
              << "; bh9.seq_length(1): " << bh9.seq_length(1) << "\n";
    std::string ss = bh9.seq_name(0);
    std::cout << "ss = bh9.seq_name(0): " << ss << "\n";
    std::cout << "bh9.name2id(CHROMOSOME_III): " << bh9.name2id("CHROMOSOME_III") << "\n";
    std::cout << "bh9 Header: \n" << bh9.header_txt() << "\n";
    int d1 = 5 / 3;
    std::cout << "HHHH: " << d1 << " -- " << 5%3 << "\n";
    
    std::cout << "bh9.get_sample_name: " << bh9.get_sample_name() << "\n";
    // std::cout << "bh7.get_sample_name: " << bh7.get_sample_name() << "\n";  // no RG in header, throw error here

    sam_close(fp);
    return 0;
}

