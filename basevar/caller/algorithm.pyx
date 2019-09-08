"""
This module contain some main algorithms of BaseVar
"""
from scipy.stats.distributions import norm

from basevar.io.htslibWrapper cimport kt_fisher_exact

cdef extern from "math.h":
    double log10(double)


cdef void EM(double* init_allele_freq,
             double* ind_allele_likelihood,
             double* marginal_likelihood,
             double* expect_allele_prob,
             int nsample,
             int ntype,
             int iter_num,
             double epsilon):

    em(init_allele_freq, ind_allele_likelihood, marginal_likelihood,
       expect_allele_prob, nsample, ntype, iter_num, epsilon)

    return


cdef double ref_vs_alt_ranksumtest(bytes ref_base, list alt_base, char **bases, int *info, int data_size):
    """Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.

    ``bases`` : A string array
             A tuple content pair-data for sample_base with other.
             
    ``info`` : A integer array
        ``info`` is the same size as ``bases`` for record information for bases
        e.g: ``bases`` = sample_base and ``info`` is mapqs (mapping quality)

    """
    cdef double* ref = <double*>(malloc(data_size * sizeof(double)))
    cdef double* alt = <double*>(malloc(data_size * sizeof(double)))

    cdef int i = 0
    cdef int size_ref = 0
    cdef int size_alt = 0
    cdef char *b
    cdef int d
    for i in range(data_size):

        b, d = bases[i], info[i]
        if b[0] in ['N', '-', '+']:
            continue

        if b == ref_base:
            ref[size_ref] = d
            size_ref += 1

        elif b in alt_base:
            alt[size_alt] = d
            size_alt += 1

    if size_ref == 0 or size_alt == 0:
        free(ref)
        free(alt)
        # -1 represent to None
        return -1.0

    cdef double* x = <double*>(malloc(size_ref * sizeof(double)))
    cdef double* y = <double*>(malloc(size_alt * sizeof(double)))
    memcpy(x, ref, size_ref * sizeof(double))
    memcpy(y, alt, size_alt * sizeof(double))

    free(ref)
    free(alt)

    cdef double z = RankSumTest(x, size_ref, y, size_alt)
    cdef double pvalue = 2 * norm.sf(abs(z))

    cdef double phred_scale_value
    if pvalue == 1.0:
        phred_scale_value = 0.0
    elif pvalue > 0:
        phred_scale_value = -10 * log10(pvalue)
    else:
        phred_scale_value = 10000.0

    free(x)
    free(y)

    return phred_scale_value


cdef tuple strand_bias(bytes ref_base, list alt_bases, char **bases, char *strands, int size):
    """
    A method for calculating the strand bias of REF_BASE and ALT_BASE

    :param ref_base:  char, required
        The reference base

    :param alt_bases: array like, required
        A list of alt bases

    :param sample_base: array-like, required
        A list of bases cover this position

    :param strands: array-like, equired
        '+' or '-' strand for each base in ``sample_base``

    :return: list-like
        FS, ref_fwd, ref_rev, alt_fwd, alt_rev
    """
    cdef int ref_fwd = 0, ref_rev = 0, alt_fwd = 0, alt_rev = 0

    cdef int i = 0
    # for s, b in zip(strands, bases):
    for i in range(size):

        # ignore "N" or indels
        if bases[i][0] in ['N', '-', '+']:
            # `bases[i]` is char*
            continue

        if strands[i] == '+':
            if bases[i] == ref_base:
                ref_fwd += 1

            elif bases[i] in alt_bases:
                alt_fwd += 1

        elif strands[i] == '-':
            if bases[i] == ref_base:
                ref_rev += 1
            elif bases[i] in alt_bases:
                alt_rev += 1

        else:
            raise ValueError('[ERROR] Get strange strand symbol: "%s"' % chr(strands[i]))

    cdef double left_p, right_p, twoside_p, fs

    # exact_fisher_test from htslib
    kt_fisher_exact(ref_fwd, ref_rev, alt_fwd, alt_rev, &left_p, &right_p, &twoside_p)
    if twoside_p == 1.0:
        fs = 0.0
    elif twoside_p > 0.0:
        fs = -10 * log10(twoside_p)
    else:
        fs = 10000.0

    # Strand bias estimated by the Symmetric Odds Ratio test
    # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    cdef double sor = float(ref_fwd * alt_rev) / (ref_rev * alt_fwd) if ref_rev * alt_fwd > 0 else 10000.0
    return (fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev)
