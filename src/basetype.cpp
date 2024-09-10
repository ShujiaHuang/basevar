/**
 * @file basetype.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */

#include <cctype>   // use toupper()

#include <htslib/kfunc.h>

#include "basetype.h"
#include "algorithm.h"
#include "external/combinations.h"

#include "utils.h"  // join(), sum()

//////////////////////////////////////////////////////////////
//// The codes for the member function of BaseType class /////
//////////////////////////////////////////////////////////////
BaseType::BaseType(const BatchInfo *smp_bi, double min_af) : _only_call_snp(true) {

    // set common values
    for (size_t i(0); i < BASES.size(); ++i) {
        _B_IDX[BASES[i]] = i;  // set [A, C, G, T] char map to array index
        _depth[BASES[i]] = 0;  // inital the base depth to 0.
    }

    _min_af      = min_af;
    _total_depth = 0;

    _ref_id   = smp_bi->ref_id;
    _ref_pos  = smp_bi->ref_pos;
    _ref_base = smp_bi->ref_base;
    _alt_bases.clear();

    // Initialized the array to {0, 0, 0, 0}, which set allele likelihood for [A, C, G, T]
    std::vector<double> allele_lh(BASES.size(), 0);
    _ind_allele_likelihood.reserve(smp_bi->n);
    _qual_pvalue.reserve(smp_bi->n);

    double epsilon; // base error probability
    char fb;
    for (size_t i(0); i < smp_bi->n; ++i) {

        epsilon = exp((smp_bi->align_base_quals[i] - 33) * MLN10TO10);
        _qual_pvalue.push_back(1.0 - epsilon);

        fb = smp_bi->align_bases[i][0];  // a string, get first base
        if (fb != 'N' && ((fb != '+' && fb != '-') || !_only_call_snp)) { 
            // ignore all the 'N' bases and indels if only call snp

            if (smp_bi->align_bases[i].size() != 1)
                throw std::runtime_error("[ERROR] Why dose the size of aligned base is not 1? Check: " + 
                                         smp_bi->align_bases[i]); 
            
            _depth[fb]++;  // record the depth for read base: [A, C, G, T]
            _total_depth++;

            for(size_t j(0); j < BASES.size(); ++j) {
                // convert the quality phred scale to be the base confident probabilty value
                allele_lh[j] = (fb == BASES[j]) ? 1.0 - epsilon : epsilon / 3;
            }

            // A 2d-array, n x 4 matrix, n is the non-N and non-Indels' sample size. (n maybe not equality to smp_bi->n)
            // 也就是说这个 likelihood 数组将只保留有覆盖的位点信息，避免无效计算。因此 _ind_allele_likelihood 的长度较小
            // 但无所谓，因为除了该值不会传到 basetype classe 之外，仅仅只用于计算突变 
            _ind_allele_likelihood.push_back(allele_lh);
        }
    }
}

BaseType::BaseType(const BaseType &b) {

    this->_only_call_snp = b._only_call_snp;
    this->_ref_id        = b._ref_id;
    this->_ref_pos       = b._ref_pos;
    this->_ref_base      = b._ref_base;
    this->_alt_bases     = b._alt_bases;

    this->_min_af      = b._min_af;
    this->_var_qual    = b._var_qual;
    this->_af_by_lrt   = b._af_by_lrt;
    this->_qual_pvalue = b._qual_pvalue;
    this->_ind_allele_likelihood = b._ind_allele_likelihood;

    this->_depth       = b._depth;
    this->_total_depth = b._total_depth;
    this->_B_IDX       = b._B_IDX;
}

// We do not have free any data before assige, no need this function.
// BaseType &BaseType::operator=(const BaseType &b) { 

//     this->_only_call_snp = b._only_call_snp;
//     this->_ref_id        = b._ref_id;
//     this->_ref_pos       = b._ref_pos;
//     this->_ref_base      = b._ref_base;
//     this->_alt_bases     = b._alt_bases;

//     this->_min_af      = b._min_af;
//     this->_var_qual    = b._var_qual;
//     this->_af_by_lrt   = b._af_by_lrt;
//     this->_qual_pvalue = b._qual_pvalue;
//     this->_ind_allele_likelihood = b._ind_allele_likelihood;

//     this->_depth       = b._depth;
//     this->_total_depth = b._total_depth;
//     this->_B_IDX       = b._B_IDX;

//     return *this;
// }

std::vector<double> BaseType::_set_allele_initial_freq(const std::vector<char> &bases) {
    // bases 数组中 A,C,G,T 这四个碱基最多只能各出现一次
    // Initialized the array to {0, 0, 0, 0}, which set initial observed allele likelihood for [A, C, G, T]
    std::vector<double> obs_allele_freq(BASES.size(), 0);
    if (this->_total_depth > 0) {
        // computed by base count
        for (auto b: bases) obs_allele_freq[_B_IDX[b]] = this->_depth[b] / (double)(this->_total_depth);
    }

    return obs_allele_freq;  // 1 x 4 vector. The allele frequence of [A, C, G, T]
}

/**
 * @brief Calculate population likelihood for all the combination of bases
 * 
 * @param bases A 1-d array. An array subset of bases from [A, C, G, T] 
 * @param n     The combination number. n must less or equal to the length of ``bases``
 * 
 * @return AA   AA.bc: An array of combinamtion bases
 *              AA.lr: Likelihood of ``bc``
 * 
 */
AA BaseType::_f(const std::vector<char> &bases, int n) {

    AA data;
    Combinations<char> c(bases, n);
    std::vector<std::vector<char>> cbs_v = c.get();  // combination bases vector
    for (size_t i = 0; i < cbs_v.size(); i++) {      // 循环该位点每一种可能的碱基组合

        std::vector<double> obs_allele_freq = this->_set_allele_initial_freq(cbs_v[i]);
        if (ngslib::sum(obs_allele_freq) == 0) // Empty coverage for this type of combination, skip.
            throw std::runtime_error("The sum of frequence of active bases must always > 0. Check: " + 
                                     ngslib::join(cbs_v[i], ",") + " - " + ngslib::join(obs_allele_freq, ","));
        
        std::vector<double> log_marginal_likelihood;
        // The value of 'obs_allele_freq' and 'log_marginal_likelihood' will be updated in EM process.
        EM(_ind_allele_likelihood, obs_allele_freq, log_marginal_likelihood);
        double sum_log_marginal_likelihood = ngslib::sum(log_marginal_likelihood);

        data.bc.push_back(cbs_v[i]);
        data.bp.push_back(obs_allele_freq);
        data.lr.push_back(sum_log_marginal_likelihood);
    }
    
    return data;
}

/**
 * @brief The main function for likelihood ratio test
 * 
 * @param specific_bases Calculating LRT for specific base combination
 * 
 */
void BaseType::lrt(const std::vector<char> &specific_bases) {

    if (_total_depth == 0) return;
    
    std::vector<char> active_bases;
    for (auto b: specific_bases) {
        // Get active bases which count frequence > _min_af
        if (_depth[b] / _total_depth >= _min_af)
            active_bases.push_back(b);
    }

    if (active_bases.size() == 0) return;

    // init. Base combination of active_bases
    AA var = _f(active_bases, active_bases.size());  // F4

    double chi_sqrt_value = 0;
    std::vector<double> active_bases_freq = var.bp[0];
    double lr_alt = var.lr[0];  // f4

    // Find candinate altnative alleles
    for (size_t n = active_bases.size() - 1; n > 0; --n) {
        var = _f(active_bases, n);
        std::vector<double> lrt_chivalue;
        for (size_t j(0); j < var.lr.size(); ++j) {
            lrt_chivalue.push_back(2 * (lr_alt - var.lr[j]));
        }
        size_t i_min = ngslib::argmin(lrt_chivalue.begin(), lrt_chivalue.end());

        lr_alt = var.lr[i_min];
        chi_sqrt_value = lrt_chivalue[i_min];
        if (chi_sqrt_value < LRT_THRESHOLD) {
            // Take the null hypothesis and continue
            active_bases = var.bc[i_min];
            active_bases_freq = var.bp[i_min];
        } else {
            // Take the alternate hypothesis
            break;
        }
    }

    char upper_ref_base = toupper(this->_ref_base[0]);  // Only call SNP: only get the first base for SNP
    for (auto b: active_bases) {
        if (b != upper_ref_base) {
            this->_alt_bases.push_back(b);
            this->_af_by_lrt[b] = active_bases_freq[_B_IDX[b]];
        }
    }

    // Todo: improve the calculation method for var_qual
    if (!this->_alt_bases.empty()) {

        double r = this->_depth[active_bases[0]] / (double)(this->_total_depth);
        if ((active_bases.size() == 1) && (this->_total_depth > 10) && (r > 0.5)) {
            // mono-allelelic
            this->_var_qual = 5000.0;
        } else {
            // 'chi2_test' may return nan, which is caused by 'chi_sqrt_value' <= 0 and means p value is 1.0.
            double chi_prob = chi2_test(chi_sqrt_value, 1);  // Python: chi_prob = chi2.sf(chi_sqrt_value, 1)
            if (std::isnan(chi_prob)) 
                chi_prob = 1.0;

            this->_var_qual = (chi_prob) ? -10 * log10(chi_prob) : 10000.0;
            // _var_qual will been setted as -0.0 instand of 0.0 if it's 0, because of the phred-scale formular
            if (this->_var_qual == -0.0) this->_var_qual = 0.0;
        }
    }

    return;
}

/**
 * @brief Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.
 * 
 * @param ref_base A reference base
 * @param alt_bases_string 
 * @param bases 
 * @param values 
 * @return double   Phred-scale value of ranksum-test pvalue
 * 
 * Note: There's some difference between scipy.stats.ranksums with R's wilcox.test:
 *       https://stackoverflow.com/questions/12797658/pythons-scipy-stats-ranksums-vs-rs-wilcox-test
 * 
 */
double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<int> &values)  // values 和 bases 的值是配对的，一一对应 
{
    std::vector<double> ref, alt;
    ref.reserve(values.size());  // change capacity and save time
    alt.reserve(values.size());  // change capacity and save time

    for (size_t i = 0; i < bases.size(); i++) {
        if (bases[i] == 'N' || bases[i] == '-' || bases[i] == '+') 
            continue;

        if (bases[i] == ref_base) {
            ref.push_back(values[i]);
        } else if (alt_bases_string.find(bases[i]) != std::string::npos) {
            alt.push_back(values[i]);
        }
    }
    
    double p_phred_scale_value;
    if (ref.size() > 0 && alt.size() > 0) { // not empty
        double p_value = wilcoxon_ranksum_test(ref, alt);
        p_phred_scale_value = -10 * log10(p_value);
        if (std::isinf(p_phred_scale_value)) {
            p_phred_scale_value = 10000;
        }
    } else {
        // set a big enough phred-scale value if only get REF or ALT base on this postion.
        p_phred_scale_value = 10000;
    }
    return p_phred_scale_value;
}

double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<char> &values) 
{
    std::vector<int> v(values.begin(), values.end()); 
    return ref_vs_alt_ranksumtest(ref_base, alt_bases_string, bases, v);
}

StrandBiasInfo strand_bias(const char ref_base, 
                           const std::string alt_bases_string,
                           const std::vector<char> &bases,
                           const std::vector<char> &strands)  // strands 和 bases 是配对的，一一对应 
{
    int ref_fwd = 0, ref_rev = 0;
    int alt_fwd = 0, alt_rev = 0;

    for (size_t i(0); i < bases.size(); ++i) {
        // ignore non bases or indels
        if (bases[i] == 'N' || bases[i] == '-' || bases[i] == '+') 
            continue;
        
        if (strands[i] == '+') {
            if (bases[i] == ref_base) {
                ++ref_fwd;
            } else if (alt_bases_string.find(bases[i]) != std::string::npos) {
                ++alt_fwd;
            }

        } else if (strands[i] == '-') {
            if (bases[i] == ref_base) {
                ++ref_rev;
            } else if (alt_bases_string.find(bases[i]) != std::string::npos) {
                ++alt_rev;
            }

        } else {
            throw std::runtime_error("[ERROR] Get strange strand symbol: " + ngslib::tostring(strands[i]));
        }
    }

    // 如果是'全 Ref' 或者是 '全 ALT' 怎么办？
    double fs = -10 * log10(fisher_exact_test(ref_fwd, ref_rev, alt_fwd, alt_rev));
    if (std::isinf(fs)) {
        fs = 10000;
    } else if (fs == 0) {
        fs = 0.0;
    }

    // Strand bias estimated by the Symmetric Odds Ratio test
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    double sor = (ref_rev * alt_fwd > 0) ? (double)(ref_fwd * alt_rev) / (double)(ref_rev * alt_fwd): 10000;

    StrandBiasInfo sbi;
    sbi.ref_fwd = ref_fwd; sbi.ref_rev = ref_rev; 
    sbi.alt_fwd = alt_fwd; sbi.alt_rev = alt_rev;
    sbi.fs  = fs; 
    sbi.sor = sor;

    return sbi;
}


