"""
This module contain some main algorithms of BaseVar
"""
from scipy.stats.distributions import norm
import numpy as np


def EM(double[::1] init_allele_freq, double[:,::1] ind_allele_likelihood,
       int iter_num=100, double epsilon=0.001):

    cdef int nsample = ind_allele_likelihood.shape[0]
    cdef int ntype = ind_allele_likelihood.shape[1]
    cdef double[::1] marginal_likelihood = np.zeros((nsample,), dtype=np.double)
    cdef double[::1] expect_allele_prob = np.zeros((ntype,), dtype=np.double)
    em(&init_allele_freq[0], &ind_allele_likelihood[0, 0], &marginal_likelihood[0],
       &expect_allele_prob[0], nsample, ntype, iter_num, epsilon)

    return np.asarray(marginal_likelihood), np.asarray(expect_allele_prob)


def ref_vs_alt_ranksumtest(ref_base, alt_base, data):
    """Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.

    ``data`` : A 2-D list,
             A tuple content pair-data for sample_base with other.
             e.g: zip(sample_base, mapq)

    """
    cdef list ref = []
    cdef list alt = []
    for b, d in data:

        if b == 'N' or b[0] in ['-', '+']:
            continue

        if b == ref_base:
            ref.append(d)

        elif b in alt_base:
            alt.append(d)

    cdef double[::1] x = np.asarray(ref, dtype=np.double)
    cdef double[::1] y = np.asarray(alt, dtype=np.double)
    cdef int nx = x.shape[0]
    cdef int ny = y.shape[0]

    if nx == 0 or ny == 0:
        return np.nan

    cdef double z = RankSumTest(&x[0], nx, &y[0], ny)
    cdef double pvalue = 2 * norm.sf(abs(z))
    cdef double phred_scale_value = round(-10 * np.log10(pvalue), 3)
    if phred_scale_value == np.inf:
        phred_scale_value = 10000.0

    return phred_scale_value if phred_scale_value != 0 else 0.0


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

    cdef double left_p, right_p, twoside_p
    # exact_fisher_test from htslib
    kt_fisher_exact(ref_fwd, ref_rev, alt_fwd,
                    alt_rev, & left_p, & right_p, & twoside_p)
    fs = round(-10 * np.log10(twoside_p), 3)

    if fs == np.inf:
        fs = 10000.0

    if fs == 0:
        # ``fs`` will been setted as -0.0 instand of 0.0 if it's 0,
        # and I don't know why it so weird!
        fs = 0.0

    # Strand bias estimated by the Symmetric Odds Ratio test
    # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    sor = round(float(ref_fwd * alt_rev) / (ref_rev * alt_fwd), 3) \
        if ref_rev * alt_fwd > 0 else 10000.0

    return fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev
