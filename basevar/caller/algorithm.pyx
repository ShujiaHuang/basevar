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


cdef double ref_vs_alt_ranksumtest(bytes ref_base, list alt_base, list data):
    """Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.

    ``data`` : A 2-D list,
             A tuple content pair-data for sample_base with other.
             e.g: zip(sample_base, mapq)

    """
    cdef int data_size = len(data)
    cdef double* ref = <double*>(malloc(data_size * sizeof(double)))
    cdef double* alt = <double*>(malloc(data_size * sizeof(double)))

    cdef int i = 0
    cdef int size_ref = 0
    cdef int size_alt = 0
    for i in range(data_size):

        b, d = data[i]
        if b == 'N' or b[0] == '-' or b[0] == '+':
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


def strand_bias(ref_base, alt_base, bases, strands):
    """
    A method for calculating the strand bias of REF_BASE and ALT_BASE

    :param ref_base:  char, required
        The reference base

    :param alt_base: array like, required
        A list of alt bases.

    :param sample_base: array-like, required
        A list of bases cover this position

    :param strands: array-like, equired
        '+' or '-' strand for each base in ``sample_base``

    :return: list-like
        FS, ref_fwd, ref_rev, alt_fwd, alt_rev
    """
    cdef int ref_fwd = 0, ref_rev = 0, alt_fwd = 0, alt_rev = 0

    for s, b in zip(strands, bases):

        # ignore non bases or indels
        if b == 'N' or b[0] in ['-', '+']:
            continue

        if s == '+':
            if b == ref_base:
                ref_fwd += 1
            elif b in alt_base:
                alt_fwd += 1

        elif s == '-':
            if b == ref_base:
                ref_rev += 1
            elif b in alt_base:
                alt_rev += 1

        else:
            raise ValueError('[ERROR] Get strange strand symbol %s' % s)

    cdef double left_p, right_p, twoside_p, fs

    # exact_fisher_test from htslib
    kt_fisher_exact(ref_fwd, ref_rev, alt_fwd,
                    alt_rev, &left_p, &right_p, &twoside_p)

    if twoside_p == 1.0:
        fs = 0.0
    elif twoside_p > 0.0:
        fs = -10 * log10(twoside_p)
    else:
        fs = 10000.0

    # Strand bias estimated by the Symmetric Odds Ratio test
    # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    sor = float(ref_fwd * alt_rev) / (ref_rev * alt_fwd) if ref_rev * alt_fwd > 0 else 10000.0
    return fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev
