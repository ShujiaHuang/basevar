/*
 *  The main program of BaseVar.
 *
 *  Created on: Jul 30, 2018
 *      Author: Shujia Huang
 * 
 * 
 */
#include <iostream>
#include <getopt.h>
#include <ctime>

// inluding my own lib
#include "basetype_caller.h"

static const std::string AUTHOR = "Shujia Huang (hshujia@qq.com)";
static const std::string VERSION = "version 1.0.1";

static int usage() {
    // Usage discription
    std::cout << 
        "Author:  " + AUTHOR  + "\n"
        "Version: " + VERSION + "\n\n"
        "Usage: basevar <command> [options]\n\n" 
        "Commands:\n"
        "    basetype            Variants Caller\n"
        // "    VQSR                Variants Recalibrator\n"
        // "    coverage            Calculating coverage depth for the whole genome or given regions/positions\n"
        // "    merge               Merge bed/vcf files\n"
        // "    NearByIndel         Calculating and adding Nearby Indel density and indel type information for each variants in VCF.\n"
        // "    popmatrix           Create population matrix from bamfiles in specific positions.\n"
        << "\n" << std::endl;

    return 1;
}

int basetype(int argc, char *argv[]) {

    BaseTypeRunner bt;
    if (argc < 2) {
        std::cout << bt.usage() << "\n" << std::endl;
        exit(1);
    }

    bt.set_arguments(argc, argv);
    bt.run();

    return 0;
}

int main(int argc, char *argv[]) {
    clock_t start_time = clock();
    std::cout << "BaseVar: A software for calling variants efficiently "
              << "from low-pass whole genome sequencing data.\n\n";
    if (argc < 2) {
        return usage();
    }

    time_t now = time(0);
    std::cout << "Program start on " << ctime(&now) << "\n";

    int run_stat;
    std::string cmd(argv[1]);
    if (cmd == "basetype") {
        run_stat = basetype(argc-1, argv+1);
    } else {
        return usage();
    }

    /* Coding here */

    
    now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "\n** " + ct + ". Processing all done, "
              << (double)(clock() - start_time) / CLOCKS_PER_SEC 
              << " seconds elapsed in total. **\n" << std::endl;

    return 0;
}
