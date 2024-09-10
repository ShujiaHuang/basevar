// Author: Shujia Huang
// Date: 2021-08-25
#include <iostream>
#include <string>
#include <ctime>

#include "bam.h"
#include "bam_record.h"

void ret_br(ngslib::Bam &b) {

    // reading alignment record
    int read_count = 0;

    // fetch alignments
    ngslib::BamRecord al;
    while (b.next(al) >= 0) {  // -1 => on end of the file.
        ++read_count;
        std::cout << "* Read count: " << read_count
                  << " ; Read status: " << b.io_status() << "\n" << al << "\n";
    }
}

#include <vector>

int main() {
    clock_t start, end;  // time recorder
    
    using ngslib::Bam;
    using ngslib::BamRecord;

    start = clock();  

    const char *fn1 = "../data/range.cram";
    std::string fn2 = "../data/range.bam";
    std::string fn3 = "../data/xx_minimal.sam";

    Bam b0;
    Bam b1(fn1, "r");
    // Bam b2(fn2, "r");
    Bam b2("../data/range.bam", "r");
    const Bam *b3; b3 = &b1;
    b0 = b1;
std::cout << b0.header() << "\n";

std::vector<Bam> bb;  
bb.reserve(5); //  capacity changed
for(size_t i(0); i < 2; ++i) {
    bb.push_back(Bam(fn1, "r"));
    std::cout << i << " - " << bb[i] << "\n";
}
std::cout << "bb size: " << bb.size() << " 0 - " << bb[0] << "\n";
exit(1);

    if (b1.index_build() == 0)
        std::cout << "Successful generate BAI-format index for BAM "
                     "files [default].\n";

    // if (b1.index_build(1) == 0)
    //    std::cout << "Successful generate CSI-format index for BAM files.\n";

    // Bam b3 = b1;   // Not allow
    // b0 = b1;       // Not allow
    // b0 = fn1;      // Not allow

    std::cout << b0.header() << "\n";
    std::cout << b1.header() << "\n";
    std::cout << ">>> 0. The file name is: " << b0 << "\n";
    std::cout << ">>> 1. The file name is: " << b1 << "\n";
    std::cout << ">>> 2. The file name is: " << b2 << "\n";
    std::cout << ">>> 3. The file name is: " << b3 << "\t" << *b3 << "\n\n";

    std::cout << "\n** Loop all the data **\n";
    ret_br(b1);

    bool good;
    std::cout << "\n** Loop CHROMOSOME_IV the data **\n";
    good = b1.fetch("CHROMOSOME_IV");
    ret_br(b1);
    std::cout << "End loop status: " << good << "\n\n";

    std::cout << "\n** Loop CHROMOSOME_I the data **\n";
    good = b1.fetch("CHROMOSOME_I");
    ret_br(b1);
    std::cout << "End loop status: " << good << "\n\n";

    std::cout << "\n** Loop CHROMOSOME_I:1-10 the data **\n";
    good = b1.fetch("CHROMOSOME_I:1-10");
    ret_br(b1);
    std::cout << "End loop status: " << good << "\n\n";

    std::cout << "\n** Loop CHROMOSOME_I:914- the data **\n";
    good = b1.fetch("CHROMOSOME_I:914-");
    ret_br(b1);
    std::cout << "End loop status: " << good << "\n\n";

    std::cout << "\n** Loop CHROMOSOME_I:-934 the data **\n";
    good = b1.fetch("CHROMOSOME_I:-934");
    ret_br(b1);
    std::cout << "End loop status: " << good << "\n\n";

    std::cout << "\n** Loop CHROMOSOME_I:914-914 the data **\n";
    good = b1.fetch("CHROMOSOME_I:914-914");
    ret_br(b1);
    std::cout << "End loop status: " << good << "\n\n";

    time_t now = time(0);
    std::cerr << ctime(&now) << "** Process done, " << (double)(clock() - start) / CLOCKS_PER_SEC
              << " seconds elapsed **\n";

    return 0;
}