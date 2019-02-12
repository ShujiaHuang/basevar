"""
This module contain some main algorithms of BaseVar
"""
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import ranksums


def EM(allele_frequence, ind_allele_likelihood, iter_num=100, epsilon=0.001):
    """
    EM algorithm

    Parameters
    ----------

    ``allele_frequence``: 2d-array like the same as ind_allele_likelihood

    ``ind_allele_likelihood`` : 2d-array like, n x 4 matrix.

    ``iter_num`` : integer, optional
        The lager EM iteration times. default: 100

    ``epsilon`` : float, optional
        The threshold of likelihood different for EM process. default 0.001
    """
    # ind_allele_likelihood is a `n x 4` matrix. n is sample size
    ind_allele_prob, marginal_likelihood = e_step(allele_frequence, ind_allele_likelihood)
    log_marginal_likelihood = np.log(marginal_likelihood)

    # m_step
    expect_allele_freq = m_step(ind_allele_prob)

    for i in range(iter_num):

        # e_step
        ind_allele_prob, marginal_likelihood = e_step(
            np.tile(expect_allele_freq, (ind_allele_likelihood.shape[0], 1)),
            ind_allele_likelihood
        )

        # m_step
        expect_allele_freq = m_step(ind_allele_prob)

        new_log_marginal_likelihood = np.log(marginal_likelihood)
        delta = np.abs(new_log_marginal_likelihood - log_marginal_likelihood).sum()

        # Todo: be careful here!!!
        if delta < epsilon:
            break

        log_marginal_likelihood = new_log_marginal_likelihood

    # update the lastest expect_allele_freq
    expect_allele_freq = m_step(ind_allele_prob)
    return ind_allele_prob, marginal_likelihood, expect_allele_freq


def e_step(allele_freq, ind_allele_likelihood):
    """
    Calculate the posterior probability of individual allele at each site as
    the four A/C/G/T bases.
    """
    likelihood = allele_freq * ind_allele_likelihood

    # It's a 1-d array one sample per element(value)
    marginal_likelihood = likelihood.sum(axis=1)

    # Posterior probability of each individual for A/C/G/T
    ind_allele_prob = likelihood / np.tile(
        marginal_likelihood.reshape(marginal_likelihood.shape[0], 1),  # convert row => col
        (1, likelihood.shape[1])  # coloum as the same as ``likelihood``
    )

    return ind_allele_prob, marginal_likelihood


def m_step(ind_base_prob):
    """
    Update allele frequency.

    Parameters
    ----------

    ``ind_base_prob``: 2d-array like, n x 4 matrix.

    """
    sample_size = ind_base_prob.shape[0]
    # ``expect_allele_prob`` is the update of ``ind_base_prob``
    expect_allele_prob = ind_base_prob.sum(axis=0) / sample_size

    return expect_allele_prob


def ref_vs_alt_ranksumtest(ref_base, alt_base, data):
    """Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.

    ``data`` : A 2-D list,
             A tuple content pair-data for sample_base with other.
             e.g: zip(sample_base, mapq)

    ``Note`` : There's some difference between scipy.stats.ranksums
               with R's wilcox.test
               https://stackoverflow.com/questions/12797658/pythons-scipy-stats-ranksums-vs-rs-wilcox-test
    """
    ref, alt = [], []
    for b, d in data:

        if b == 'N' or b[0] in ['-', '+']:
            continue

        if b == ref_base:
            ref.append(d)

        elif b in alt_base:
            alt.append(d)

    _, pvalue = ranksums(ref, alt)
    phred_scale_value = round(-10 * np.log10(pvalue), 3)
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
    ref_fwd, ref_rev = 0, 0
    alt_fwd, alt_rev = 0, 0

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

    # Use R to calculate the strand bias by fisher exact test instand of scipy
    # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
    # m = R['matrix'](robjects.IntVector([ref_fwd, alt_fwd,
    #                                     ref_rev, alt_rev]), nrow=2)
    # pvalue = R['fisher.test'](m)[0][0]
    # fs = round(-10 * np.log10(pvalue), 3)

    # May be too slow?
    fs = round(-10 * np.log10(fisher_exact([[ref_fwd, ref_rev],[alt_fwd, alt_rev]])[1]), 3)

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
