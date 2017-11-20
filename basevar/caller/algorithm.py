"""
This module contain some main algorithms of BaseVar
"""
import sys

import numpy as np
# from scipy.stats import fisher_exact
from rpy2 import robjects
from rpy2 import rinterface

R = robjects.r


def EM(prior_prob, ind_base_likelihood, iter_num=100, epsilon=0.001):
    """
    EM Process

    Parameters
    ----------

    ``prior_prob``: array-like

    ``ind_base_likelihood`` : 2d-array like, n x 4 matrix.

    ``iter_num`` : integer, optional
        The lager EM iteration times. default: 100

    ``epsilon`` : float, optinal
        The threshold of likelihood different for EM process. default 0.001
    """
    # ind_base_likelihood is a `n x 4` matrix. n is sample size
    ind_base_likelihood, pop_likelihood = m_step(prior_prob,
                                                 ind_base_likelihood)
    log_pop_likelihood = np.log(pop_likelihood)
    for i in xrange(iter_num):

        # e_step
        base_expect_prob = e_step(ind_base_likelihood)

        # m_step
        ind_base_likelihood, pop_likelihood = m_step(
            prior_prob,
            np.tile(base_expect_prob, (prior_prob.shape[0], 1)))

        new_log_pop_likelihood = np.log(pop_likelihood)
        delta = np.abs(new_log_pop_likelihood - log_pop_likelihood).sum()
        if delta < epsilon:
            break

        log_pop_likelihood = new_log_pop_likelihood

    base_expect_prob = e_step(ind_base_likelihood)  # update the lastest base_expect_prob
    return ind_base_likelihood, pop_likelihood, base_expect_prob


def e_step(ind_base_likelihood):
    """
    Expection step: update ``ind_base_likelihood``

    Calculate the base posterior probability.

    Parameters
    ----------

    ``ind_base_likelihood``: 2d-array like, n x 4 matrix.

    """
    sample_size = ind_base_likelihood.shape[0]
    # ``base_expect_prob`` is the update of ``ind_base_likelihood``
    base_expect_prob = ind_base_likelihood.sum(axis=0) / sample_size

    return base_expect_prob


def m_step(prior_prob, ind_base_likelihood):
    """
    Maximization the likelihood step
    """
    likelihood = prior_prob * ind_base_likelihood

    # It's a 1-d array one sample per element
    pop_likelihood = likelihood.sum(axis=1)

    # Maximize the base likelihood
    ind_base_likelihood = likelihood / np.tile(
        pop_likelihood.reshape(pop_likelihood.shape[0], 1),  # convert row => col
        (1, likelihood.shape[1]))  # coloum as the same as ``likelihood``

    return ind_base_likelihood, pop_likelihood


def ref_vs_alt_ranksumtest(ref_base, alt_base, data):
    """Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.

    ``data`` : A 2-D list,
             A tuple content pair-data for sample_base with other.
             e.g: zip(sample_base, mapq)
    """
    ref, alt = [], []
    for b, d in data:

        if b == 'N':
            continue

        if b == ref_base:
            ref.append(d)

        elif b in alt_base:
            alt.append(d)

    ref = robjects.FloatVector(ref)
    alt = robjects.FloatVector(alt)

    try:
        pvalue = R['wilcox.test'](ref, alt)[2][0]
        phred_scale_value = round(-10 * np.log10(pvalue), 3)

    except rinterface.RRuntimeError:
        sys.stderr.write('[WARNING] The array number is too samll for '
                         'wilcox.test and set phred_scale_value = 0 :\n'
                         '%s\n%s\n' % (ref, alt))
        phred_scale_value = 0

    return phred_scale_value


def strand_bias(ref_base, alt_base, sample_base, strands):
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
    ref_fwd, ref_rev = 0, 0
    alt_fwd, alt_rev = 0, 0

    for s, b in zip(strands, sample_base):

        # For strand bias
        if b == 'N':
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

    # Use R to calculate the strand bias by fisher exact test instand of scipy
    # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
    m = R['matrix'](robjects.IntVector([ref_fwd, alt_fwd,
                                        ref_rev, alt_rev]), nrow=2)
    pvalue = R['fisher.test'](m)[0][0]
    fs = round(-10 * np.log10(pvalue), 3)

    # May be too slow!
    # fs = round(-10 * np.log10(fisher_exact([[ref_fwd, ref_rev],
    #                                         [alt_fwd, alt_rev]])[1]), 3)

    if fs == np.inf:
        fs = 10000.0

    # Strand bias estimated by the Symmetric Odds Ratio test
    # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    sor = round(float(ref_fwd * alt_rev) / (ref_rev * alt_fwd), 3) \
        if ref_rev * alt_fwd > 0 else 10000.0

    return fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev
