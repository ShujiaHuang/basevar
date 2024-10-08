/**
 * @file basetype_utils.h
 * 
 * @brief define the global variable and util codes specific for basetype.h
 *  
 * @author Shujia Huang
 * @date 2018-08-29
 * 
 */

#ifndef __INCLUDE_BASETYPE_UTILS_H__
#define __INCLUDE_BASETYPE_UTILS_H__

#include <getopt.h>
#include <string>
#include <vector>

static const std::string __BASETYPE_USAGE = 
    "About: Call variants and estimate allele frequency by BaseVar.\n" 
    "Usage: basevar basetype [options] <-R Fasta> <--output-vcf> <--output-cvg> [-I input] ...\n\n" 
    "optional arguments:\n" 
    "  -I, --input=FILE             BAM/CRAM file containing reads.\n"
    "  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.\n"
    "  -R, --reference FILE         Input reference fasta file.\n\n"

    "  -m, --min-af=float           Setting prior precision of MAF and skip ineffective caller positions,\n"
    "                               a typical approach involves setting it to min(0.001, 100/x), where x \n"
    "                               represents the number of input BAM files [min(0.001, 100/x)]. In most\n"
    "                               cases, users need not be overly concerned about this parameter, as it \n"
    "                               is generally handled automatically by the program.\n"
    "  -q, --mapq=INT               Only include reads with mapping quality >= INT. [10]\n"
    "  -B, --batch-count=INT        INT simples per batchfile. [200]\n" 
    "  -t, --thread=INT             Number of threads. [4]\n\n"

    "  -G, --pop-group=FILE         Calculating the allele frequency for specific population.\n" 
    "  -r, --regions=chr:start-end  Skip positions which not in these regions. This parameter could be a list\n"
    "                               of comma deleimited genome regions(e.g.: chr:start-end) or a file contain\n"
    "                               the list of regions.\n"
    "  --output-vcf FILE            Output VCF file.\n"
    "  --output-cvg FILE            Output position coverage file.\n\n"

    "  --filename-has-samplename    If the name of bamfile is something like 'SampleID.xxxx.bam', set this\n"
    "                               argrument could save a lot of time during get the sample id from BAMfile.\n"
    "  --smart-rerun                Rerun process by checking batchfiles.\n"
    "  -h, --help                   Show this help message and exit."; 

static const struct option BASETYPE_CMDLINE_LOPTS[] = {
    // Optional arguments to long style command line parameters require 'equals sign' (=). 
    // https://stackoverflow.com/questions/1052746/getopt-does-not-parse-optional-arguments-to-parameters
    {"input",           optional_argument, NULL, 'I'},
    {"align-file-list", optional_argument, NULL, 'L'},
    {"reference",       required_argument, NULL, 'R'},

    {"min-af",      optional_argument, NULL, 'm'},
    {"mapq",        optional_argument, NULL, 'q'},
    {"batch-count", optional_argument, NULL, 'B'},
    {"thread",      optional_argument, NULL, 't'},

    {"regions",     optional_argument, NULL, 'r'},
    {"positions",   optional_argument, NULL, 'p'},
    {"pop-group",   optional_argument, NULL, 'G'},  // Special parameter for calculating specific population allele frequence
    {"output-vcf",  required_argument, NULL, '1'},
    {"output-cvg",  required_argument, NULL, '2'},

    // {"output-batch-file", required_argument, NULL, 0},  not use?
    {"filename-has-samplename", no_argument, NULL, '3'},
    {"smart-rerun",             no_argument, NULL, '4'},
    {"help",                    no_argument, NULL, 'h'},

    // must set this value
    {0, 0, 0, 0}
};

struct BaseTypeARGS {
    /* Variables for all the commandline options of BaseType */
    std::vector<std::string> input_bf;  // BAM/SAM/CRAM file, a vector
    std::string in_bamfilelist;         // BAM/CRAM files list, one file per row
    std::string reference;              // Input reference fasta file

    float min_af;                       // Setting prior precision of MAF and skip uneffective caller positions
    int mapq;                           // mapping quality
    int batchcount;                     // INT simples per batchfile
    int thread_num;                     // number of threads

    std::string regions;                // Interval regions
    std::string pop_group_file;         // Specific population
    std::string output_vcf;             // Output VCF file
    std::string output_cvg;             // Output coverage file

    bool filename_has_samplename;       // sample name in file name
    bool smart_rerun;                   // Smart rerun by checking batchfiles

    // Set default argument
    BaseTypeARGS(): min_af(0.01), mapq(10), batchcount(200), thread_num(4), 
                    smart_rerun(false), filename_has_samplename(false) {}
};

// Getting the first column from input file, this it's used for getting 
// filename from input filelist.
std::vector<std::string> get_firstcolumn_from_file(const std::string fn);

/// Header for VCF
std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &addition_info, 
                              const std::vector<std::string> &samples);
std::string cvg_header_define(const std::vector<std::string> &group_info, const std::vector<char> &BASES);
void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header="#", bool is_remove_tempfile=false);

#endif