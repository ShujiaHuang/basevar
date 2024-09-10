#include <iostream>
#include <string>

#include "fasta.h"

int main() {
    using ngslib::Fasta;

    Fasta ff = "../data/ce.fa.gz";
    for (size_t i(0); i < ff.nseq(); ++i) {
        std::cout << i+1 << " - " << ff.iseq_name(i) << " - " << ff.seq_length(ff.iseq_name(i)) << "\n";
    }
    
    std::string fn = "../data/tinyfasta.fa";
    Fasta fa0;                // default constructor function
    Fasta fa(fn);             // constructor function: `Fasta(const string file_name)`
    Fasta fa1(fn.c_str()); 
    Fasta fa2 = fn;           // assignment function
    Fasta fa3 = fa;           // copy constructor
    Fasta fa_;                // default constructor function
    Fasta fa4 = fa3;

    Fasta fa5;
    fa5 = fn;
    fa5 = fa;
    fa3 = fa5;                // assignment function

    std::cout << "***** start *****\n";
    std::cout << "The path of FASTA file0: " << fa0 << std::endl;
    std::cout << "The path of FASTA file1: " << fa  << std::endl;
    std::cout << "The path of FASTA file2: " << fa2 << std::endl;
    std::cout << "The path of FASTA file3: " << fa3 << std::endl;
    std::cout << "fetch(\"ref1\", 0, 10) : " << fa.fetch("ref1", 0, 10)   << std::endl;
    std::cout << "fetch(\"ref1:1-10\"):    " << fa.fetch("ref1:1-10")     << std::endl;
    std::cout << "fa[\"ref1:1-10\"]:       " << fa["ref1:1-10"]           << std::endl;
    std::cout << "fetch(\"ref1:5\"):       " << fa.fetch("ref1:5")        << std::endl;
    std::cout << "fetch(\"ref1\", 0, 1):   " << fa1.fetch("ref1", 0, 1)   << std::endl;
    std::cout << "fa[\"ref1:1-2\"]:        " << fa1["ref1:1-2"]           << std::endl;
    std::cout << "fetch(\"ref1\", 1, 12):  " << fa2.fetch("ref1", 1, 12)  << std::endl;
    std::cout << "fetch(\"ref1\", 1, 12):  " << fa3.fetch("ref1", 1, 12)  << std::endl;
    std::cout << "fetch(\"ref1\", 1, 100): " << fa4.fetch("ref1", 1, 100) << std::endl;
    std::cout << "sequence_length(\"ref1\"): " << fa.seq_length("ref1") << std::endl;
    std::cout << "fa.has_seq(\"ref2\"): " << fa.has_seq("ref2") << std::endl;
    std::cout << "fa.has_seq(\"ref5\"): " << fa.has_seq("ref5") << std::endl;
    std::cout << "fa[\"ref2\"]:   "  << fa["ref2"] << std::endl;
    std::cout << "fa[\"ref2:1\"]: "  << fa["ref2:1"] << std::endl;
    std::cout << "fa.nseq(): "       << fa.nseq() << std::endl;
    std::cout << "fa.iseq_name(1): " << fa.iseq_name(1) << " - length: " << fa.seq_length(fa.iseq_name(1)) << std::endl;

    return 0;
}
