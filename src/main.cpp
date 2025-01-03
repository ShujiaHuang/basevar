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
#include "concat.h"

static const std::string AUTHOR = "Shujia Huang (hshujia@qq.com)";
static const std::string VERSION = "version 1.2.0";

static int usage() {
    // Usage discription
    std::cout << 
        "Program: basevar (Variant calling and allele frequency estimation from ultra low-pass WGS data)\n\n"
        "Version: " + VERSION + "\n"
        "Author:  " + AUTHOR  + "\n\n"
        "Usage: basevar <command> [options]\n\n" 
        "Commands:\n"
        "    basetype            Variants Caller\n"
        "    concat              Concatenate VCF/CVG files from the same set of samples generated by BaseVar.\n"
        // "    VQSR                Variants Recalibrator\n"
        // "    coverage            Calculating coverage depth for the whole genome or given regions/positions\n"
        // "    NearByIndel         Calculating and adding Nearby Indel density and indel type information for each variants in VCF.\n"
        // "    popmatrix           Create population matrix from bamfiles in specific positions.\n"
        << "\n" << std::endl;

    return 1;
}

// call variant
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

// Merge VCF/CVG files to a final big one
int concat(int argc, char *argv[]) {
    return concat_runner(argc, argv);
}

int main(int argc, char *argv[]) {
    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);
    
    if (argc < 2) {
        return usage();
    }

    time_t now = time(0);
    std::cout << "Program start on " << ctime(&now) << "\n";

    int run_stat;
    std::string cmd(argv[1]);
    if (cmd == "basetype") {
        run_stat = basetype(argc-1, argv+1);
    } else if (cmd == "concat") {
        run_stat = concat(argc-1, argv+1);
    } else {
        return usage();
    }

    now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "\n** " + ct + ". Processes are all done, "
              << difftime(now, real_start_time) << " (CPU time: "
              << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
              << ") seconds elapsed in total. **\n" << std::endl;

    return run_stat;
}
