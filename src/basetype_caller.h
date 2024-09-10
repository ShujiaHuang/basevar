/**
 * @file basetype_caller.h
 * 
 * @brief BaseType main caller functions
 *  
 * @author Shujia Huang
 * @date 2018-08-29
 * 
 */

#ifndef __INCLUDE_BASETYPE_CALLER_H__
#define __INCLUDE_BASETYPE_CALLER_H__

#include <map>

#include "fasta.h"
#include "bam.h"
#include "utils.h"
#include "external/robin_hood.h"

#include "basetype.h"
#include "basetype_utils.h"

static const bool IS_DELETE_CACHE_BATCHFILE = true;

// Mainly use for create batchfile
struct AlignBaseInfo {
    std::string ref_id;
    uint32_t ref_pos;
    std::string ref_base;   // reference base

    std::string read_base;  // read base
    int mapq;               // mapping quality
    int rpr;                // read position rank
    char map_strand;        // mapping reference strand, should be one of '*', '-', '+'
    char read_base_qual;    // read base quality (get mean quality of read if Indels, 
                            // I don't care about Indels for NIPT data)
};

/**
 * @brief BaseTypeRunner class
 */
class BaseTypeRunner {

private:
    BaseTypeARGS *_args;                                        // Commandline options
    std::vector<std::string> _samples_id;                       // sample ID of alignment files (BAM/CRAM/SAM)
                                                                // `_samples_id` and `input_bf` have the same order 
    std::map<std::string, std::vector<size_t>> _groups_idx;     // sample group: group => samples index
    std::vector<ngslib::GenomeRegionTuple> _calling_intervals;  // vector of calling regions

    // templary output files
    std::vector<std::string> _sub_out_vcf, _sub_out_cvg;

    void _get_bamfile_list();
    void _get_calling_interval();  // load the calling region from input
    void _get_sample_id_from_bam();
    void _get_popgroup_info();
    ngslib::GenomeRegionTuple _make_gregion_tuple(std::string gregion);

    // For variant calling
    void _variant_caller_process();
    std::vector<std::string> _create_batchfiles(const ngslib::GenomeRegionTuple &genome_region, 
                                                const std::string bf_prefix);
    void _variants_discovery(const std::vector<std::string> &batchfiles, 
                             const ngslib::GenomeRegionTuple &genome_region,
                             const std::string sub_vcf_fn,
                             const std::string sub_cvg_fn);

    BaseTypeRunner(const BaseTypeRunner &b) = delete;             // reject using copy constructor (C++11 style).
    BaseTypeRunner &operator=(const BaseTypeRunner &b) = delete;  // reject using copy/assignment operator (C++11 style).

public:
    ngslib::Fasta reference;  // public variable

    // default constructor
    BaseTypeRunner() : _args(NULL) {}
    BaseTypeRunner(int cmdline_argc, char *cmdline_argv[]) { set_arguments(cmdline_argc, cmdline_argv); }
    
    // Destroy the malloc'ed BasTypeArgs structure
    ~BaseTypeRunner(){ if(_args){delete _args; _args = NULL;} }

    // Common functions
    std::string usage() const {return __BASETYPE_USAGE;}
    void set_arguments(int cmdline_argc, char *cmdline_argv[]);
    void print_calling_interval();

    // Run the variant calling process
    void run();

};  // BaseTypeRunner class

typedef robin_hood::unordered_map<uint32_t, AlignBaseInfo> PosMap;           // give a short name to this type
typedef std::vector<PosMap> PosMapVector;

// This function is only used by BaseTypeRunner::_create_batchfiles
bool __create_a_batchfile(const std::vector<std::string> batch_align_files,  // Not a modifiable value
                          const std::vector<std::string> batch_sample_ids,   // Not a modifiable value
                          const std::string &fa_seq,                         // Not a modifiable value
                          const ngslib::GenomeRegionTuple &genome_region,    // 切割该区间
                          const int mapq_thd,                                // mapping quality threshold
                          const std::string output_batch_file);              // output batchfile name

bool __fetch_base_in_region(const std::vector<std::string> &batch_align_files,
                            const std::string &fa_seq,                   
                            const int mapq_thd,
                            const ngslib::GenomeRegionTuple target_genome_region,  // 获取该区间内的 read
                            PosMapVector &batchsamples_posinfomap_vector);         // 信息存入该变量

void __seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,  // ngslib::BamRecord include by 'bam.h'
                     const std::string &fa_seq,
                     const ngslib::GenomeRegionTuple target_genome_region,    // 获取该区间内所有位点的碱基比对信息，该参数和 '__fetch_base_in_region' 中一样 
                     PosMap &sample_posinfo_map);

void __write_record_to_batchfile(const PosMapVector &batchsamples_posinfomap_vector, 
                                 const std::string &fa_seq,
                                 const ngslib::GenomeRegionTuple target_genome_region,  // 该参数和 __seek_position 中一样 
                                 BGZF *obf);

// A unit for calling variants and let it run in a thread.
bool _variant_calling_unit(const std::vector<std::string> &batchfiles, 
                           const std::vector<std::string> &sample_ids,
                           const std::map<std::string, std::vector<size_t>> & group_smp_idx,
                           const double min_af,
                           const std::string region,  // genome region format like samtools
                           const std::string tmp_vcf_fn,
                           const std::string tmp_cvg_fn);

bool _basevar_caller(const std::vector<std::string> &smp_bf_line_vector, 
                     const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                     double min_af,
                     size_t n_sample, 
                     BGZF *vcf_hd, 
                     BGZF *cvg_hd);

const BaseType __gb(const BatchInfo *smp_bi, 
              const std::vector<size_t> group_idx, 
              const std::vector<char> &basecombination, 
              double min_af);

const BatchInfo __get_group_batchinfo(const BatchInfo *smp_bi, const std::vector<size_t> group_idx);

void _out_vcf_line(const BaseType &bt, 
                   const std::map<std::string, BaseType> &group_bt, 
                   const BatchInfo *smp_bi, 
                   BGZF *vcf_hd);

void _out_cvg_line(const BatchInfo *smp_bi, 
                   const std::map<std::string, std::vector<size_t>> & group_smp_idx, 
                   BGZF *cvg_hd);

typedef std::tuple<std::string, std::map<char, int>> IndelTuple;
IndelTuple __base_depth_and_indel(const std::vector<std::string> &align_bases);

#endif
