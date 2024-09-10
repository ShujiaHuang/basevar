/**
 * @file basetype.h
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */
#ifndef __INCLUDE_BASETYPE_H__
#define __INCLUDE_BASETYPE_H__

#include <iostream>
#include <cstdint> // uint32_t
#include <cmath>   // use exp() function
#include <string>
#include <vector>
#include <map>

static const std::vector<char> BASES = {'A', 'C', 'G', 'T'}; // 定义这个值，限定 likelihood 数组中碱基似然值的存放顺序
static const double MLN10TO10   = -0.23025850929940458;      // ln(10)/10，把 phred-value 换成 e 为底，方便调用 exp()
static const int LRT_THRESHOLD  = 24;                        // 24 corresponding to a chi-pvalue of 10^-6
static const int QUAL_THRESHOLD = 20;                        // -10 * lg(10^-2)

// Mainly use for basevar variant callling
struct BatchInfo {

    size_t n;

    std::string ref_id;
    std::string ref_base;
    uint32_t    ref_pos;
    uint32_t    depth; // The coverage depth on ref_pos

    std::vector<std::string> align_bases;
    std::vector<char>        align_base_quals;

    std::vector<int>  mapqs;
    std::vector<char> map_strands;
    std::vector<int>  base_pos_ranks;

    // Set default argument
    BatchInfo(): n(0), ref_pos(0), depth(0) {}
};

/**
 * @brief Define a structure for recording allele information return by _f() 
 * in this class.
 * 
 */
typedef struct {
    std::vector<std::vector<char>> bc;
    std::vector<std::vector<double>> bp;
    std::vector<double> lr;
} AA;

typedef struct {
    int ref_fwd, ref_rev;
    int alt_fwd, alt_rev;
    double fs;
    double sor;
} StrandBiasInfo;

// A class for calculate the base probability
class BaseType {

private:
    bool _only_call_snp;

    std::string _ref_id;
    std::string _ref_base;
    uint32_t _ref_pos;
    std::vector<char> _alt_bases; // only call SNP for now   
    double _var_qual;
    double _min_af;

    // Estimated SNP allele frequency by EM and LRT
    std::map<char, double> _af_by_lrt;

    // [A, C, G, T] likelihood vector for echo individual
    std::vector<std::vector<double>> _ind_allele_likelihood; // 2d-array, n x 4 matrix, n is sample size.
    std::vector<double> _qual_pvalue; // n x 1 matrix, n is sample size.
    std::map<char, size_t> _B_IDX;    // A map for recrod the BASE => index
    
    int _total_depth;               // sum depth of ACGT
    std::map<char, double> _depth;  // allele depth of [A, C, G, T], double 是为了方便做除法

    // init the base likelihood by input bases
    std::vector<double> _set_allele_initial_freq(const std::vector<char> &bases);

    // Calculate population likelihood for all the combination of bases
    AA _f(const std::vector<char> &bases, int n);

public:
    // Constructor
    BaseType(){};
    BaseType(const BatchInfo *smp_bi, double af);
    ~BaseType(){};

    BaseType(const BaseType &b);  // copy constructor
    // BaseType &operator=(const BaseType &b); // We do not have free any data before assige, no need this function;

    // The main function for likelihood ratio test
    /**
     * @brief 
     * 
     * @param specific_bases a 1d-array. [optional]
     * Calculating LRT for specific base combination if provided.
     */
    void lrt() { /* default */ lrt(BASES); }
    void lrt(const std::vector<char> &specific_bases);

    const bool is_only_snp() const { return this->_only_call_snp; }
    const std::string &get_ref_id() const { return this->_ref_id; };
    const uint32_t &get_ref_pos() const { return this->_ref_pos; };
    const std::string &get_ref_base() const { return this->_ref_base; };
    const std::vector<char> &get_alt_bases() const { return this->_alt_bases; };
    const std::vector<double> &get_qual_pvalue() const { return this->_qual_pvalue; }

    const double get_var_qual() const { return this->_var_qual; }
    const int get_total_depth() const { return this->_total_depth; }

    const double get_base_depth(char b) const {
        // operator[] doesn't have a 'const' qualifier in std::map. Use 'at' instead in C++11
        // https://stackoverflow.com/questions/42095642/error-passing-const-stdmapint-int-as-this-argument-discards-qualifiers
        double d;
        try {
            d = this->_depth.at(b);
        } catch (const std::out_of_range &ex) {
            std::string what_ex(ex.what());
            throw std::runtime_error("[ERROR] out_of_range:: " + what_ex + " '"+ b + "' not found.");
        }
        return d; 
    }

    const double get_lrt_af(char b) const { 
        double lrt_af;
        try {
            lrt_af = this->_af_by_lrt.at(b);
        } catch (const std::out_of_range &ex) {
            std::string what_ex(ex.what());
            throw std::runtime_error("[ERROR] out_of_range:: " + what_ex + " '"+ b + "' not found.");
        }
        return lrt_af; 
    }

}; // BaseType class

double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<int> &values);

double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<char> &values);

StrandBiasInfo strand_bias(const char ref_base, 
                           const std::string alt_bases_string,
                           const std::vector<char> &bases,
                           const std::vector<char> &strands);

#endif
