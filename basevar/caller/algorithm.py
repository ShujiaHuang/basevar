"""
This module contain some main algorithms of BaseVar
"""
import numpy as np
from scipy.stats import fisher_exact


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
    ind_base_likelihood, pop_likelihood = m_step(ind_base_likelihood)
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


def strand_bias(ref_base, alt_base, sample_base, strands):
    """
    A method for calculating the strand bias of REF_BASE and ALT_BASE

    :param ref_base:  char, required
        The reference base

    :param alt_base: array like, required
        A list of alt bases.

    :param sample_base: array-like, required
        A list of bases which cover this position

    :param strands: array-like, equired
        '+' or '-' strand for each base in ``sample_base``

    :return: list-like
        FS, ref_fwd, ref_rev, alt_fwd, alt_rev
    """
    ref_fwd, ref_rev = 0, 0
    alt_fwd, alt_rev = 0, 0
    for k, b in enumerate(sample_base):

        # For strand bias
        if b == 'N': continue
        if strands[k] == '+':
            if b == ref_base.upper():
                ref_fwd += 1
            elif b in alt_base:
                alt_fwd += 1

        elif strands[k] == '-':
            if b == ref_base.upper():
                ref_rev += 1
            elif b in alt_base:
                alt_rev += 1

    # Strand bias by fisher exact test
    # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
    fs = round(-10 * np.log10(
        fisher_exact([[ref_fwd, ref_rev], [alt_fwd, alt_rev]])[1]), 3)

    return fs, ref_fwd, ref_rev, alt_fwd, alt_rev
