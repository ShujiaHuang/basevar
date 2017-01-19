"""
This module contain functions of EM algorithm and Base genotype.

Author: Shujia Huang
Date : 2016-12-16
Update : 2017-01-03

"""
import itertools   # Use the combinations function
import numpy as np
from scipy.stats import chisqprob

from utils import CommonParameter

class BaseType(object):

    def __init__(self, ref_base, bases, quals, cmm=CommonParameter()):
        """A class for calculate the base probability

        Parameters
        ----------

        ``ref_base``: A char, required
            The reference base

        ``bases``: A array like, required
            A list of base for samples.

        ``quals``: An array like, required
            Base quality for ``bases``. The same size with ``bases``
            Cause: The ``quals`` is an integer array which has be converted
                by phred-scale

        ``cmm``: A ``CommonParameter`` class value, optional

        """
        bases = np.array(bases)
        quals = np.array(quals)

        self._alt_bases = []
        self._var_qual = 0 # init the variant quality
        self._ref_base = ref_base
        self.cmm = cmm

        #The prior probability for echo base
        self.prior_prob = []
        self.eaf = {} ## estimated allele frequency by EM
        self.depth = {b:0 for b in self.cmm.BASE}

        self.qual_pvalue = 1.0 - np.exp(self.cmm.MLN10TO10 * quals)
        for i, b in enumerate(bases):
            ## individual per row for [A, C, G, T] and the prior probability is
            if b != 'N':  # ignore all the 'N' base sample
                self.prior_prob.append([self.qual_pvalue[i]
                                        if b == t else (1.0-self.qual_pvalue[i])/3
                                        for t in self.cmm.BASE])

            # ignore all bases which not match ``cmm.BASE``
            if b in self.depth:
                self.depth[b] += 1

        self.prior_prob = np.array(self.prior_prob)
        self.total_depth = float(sum(self.depth.values()))

        return

    def ref_base(self):
        return self._ref_base

    def alt_bases(self):
        return self._alt_bases

    def var_qual(self):
        return self._var_qual

    def debug(self):
        print self.ref_base(), self.alt_bases(), self.var_qual(), self.depth, self.eaf

    def __set_base_likelihood(self, bases):
        """
        init the base likelihood by bases

        ``bases``: a list like
        """
        total_depth = float(sum([self.depth[b] for b in bases]))
        base_likelihood = np.zeros(len(self.cmm.BASE))  # [A, C, G, T] set to 0.0

        if total_depth > 0:
            for b in bases:
                base_likelihood[self.cmm.BASE2IDX[b]] = self.depth[b]/total_depth

        return np.array(base_likelihood)

    def _f(self, bases, n):
        """
        Calculate the population likelihood for all bases combination

        Parameters
        ----------

        ``bases``: 1D array like
            A list of bases in [A, C, G, T]

        ``n``: Integer
            The combination number. n must less the length of ``bases``

        Return
        ------
        ``bc`` : base array of combination
        ``lr`` : Likelihood of ``bc``

        Example
        -------

            >>> import itertools
            >>> bases = ['A', 'C', 'G', 'T']
            >>> bc=[i for i in itertools.combinations(bases,3)]
            >>> bc
            ... [('A', 'C', 'G'), ('A', 'C', 'T'), ('A', 'G', 'T'), ('C', 'G', 'T')]

        """
        bc = []
        lr = []
        bp = []

        for b in [i for i in itertools.combinations(bases, n)]:

            init_likelihood = self.__set_base_likelihood(b)

            if sum(init_likelihood) == 0: continue  ## No covert

            _, pop_likelihood, base_expected_prob = self.em(
                np.tile(init_likelihood, (self.prior_prob.shape[0],1)))

            bc.append(b)
            lr.append(np.log(pop_likelihood).sum()) # sum the marginal prob
            bp.append(base_expected_prob)

        return bc, lr, bp

    def lrt(self):
        """likelihood ratio test.
        """

        if self.total_depth == 0: return
        # init likelihood as tetra-allelic
        bc4, lr_null, bp = self._f(self.cmm.BASE, 4)
        base_frq = bp[0]
        lr_alt = lr_null[0]
        bases = self.cmm.BASE

        chi_sqrt_value = 0
        for n in (3,2,1):

            bc, lr_null, bp = self._f(bases, n)
            lrt_chivalue = 2.0 * (lr_alt - lr_null)
            i_argmin = np.argmin(lrt_chivalue)

            lr_alt = lr_null[i_argmin]
            chi_sqrt_value = lrt_chivalue[i_argmin]
            if chi_sqrt_value < self.cmm.LRT_THRESHOLD:
                # Take the null hypothesis and continue
                bases = bc[i_argmin]
                base_frq = bp[i_argmin]
            else:
                # Take the alternate hypothesis
                break

        self._alt_bases = [b for b in bases if b != self._ref_base]
        self.eaf = {b:'%f' % round(base_frq[self.cmm.BASE2IDX[b]], 6)
                    for b in bases if b != self._ref_base}

        ## calculate the variant quality
        if len(self._alt_bases):

            if len(bases) == 1 and self.depth[bases[0]] / self.total_depth > 0.5:
                # mono-allelelic
                self._var_qual = 5000

            else:
                chi_prob = chisqprob(chi_sqrt_value, 1)
                self._var_qual = round(-10 * np.log10(chi_prob), 2) if chi_prob else 10000

        return

    def em(self, ind_base_likelihood, iter_num=100, epsilon=0.001):
        """
        EM Process

        Parameters
        ----------

        ``ind_base_likelihood`` : 2d-array like, n x 4 matrix.

        ``iter_num`` : integer, optional
            The lager EM iteration times. default: 100

        ``epsilon`` : float, optinal
            The threshold of likelihood different for EM process. default 0.001
        """
        # ind_base_likelihood is a `n x 4` matrix. n is sample size
        ind_base_likelihood, pop_likelihood = self.m_step(ind_base_likelihood)
        log_pop_likelihood = np.log(pop_likelihood)
        base_expect_prob = []
        for i in xrange(iter_num):

            # e_step
            base_expect_prob = self.e_step(ind_base_likelihood)

            # m_step
            ind_base_likelihood, pop_likelihood = self.m_step(
                np.tile(base_expect_prob, (self.prior_prob.shape[0],1)))

            new_log_pop_likelihood = np.log(pop_likelihood)

            delta = np.abs(new_log_pop_likelihood - log_pop_likelihood).sum()
            if delta < epsilon:
                break

            log_pop_likelihood = new_log_pop_likelihood

        base_expect_prob = self.e_step(ind_base_likelihood)  # update the lastest base_expect_prob
        return ind_base_likelihood, pop_likelihood, base_expect_prob

    def e_step(self, ind_base_likelihood):
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

    def m_step(self, ind_base_likelihood):
        """
        Maximization the likelihood step
        """
        likelihood = self.prior_prob * ind_base_likelihood

        # It's a 1-d array one sample per element
        pop_likelihood = likelihood.sum(axis=1)

        # Maximize the base likelihood
        ind_base_likelihood = likelihood / np.tile(
            pop_likelihood.reshape(pop_likelihood.shape[0],1), # convert row => col
            (1, likelihood.shape[1])) # coloum as the same as ``likelihood``

        return ind_base_likelihood, pop_likelihood



