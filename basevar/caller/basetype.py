"""
This module contain functions of EM algorithm and Base genotype.
"""
import itertools  # Use the combinations function
import numpy as np
from scipy.stats.distributions import chi2

from .algorithm import EM


class BaseType(object):

    def __init__(self, ref_base, bases, quals, cmm=None):
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

        ``cmm``: A ``CommonParameter`` class value, required

        """

        self._alt_bases = []
        self._var_qual = 0  # init the variant quality
        self._ref_base = ref_base
        self.cmm = cmm

        # The allele likelihood for echo individual
        self.ind_allele_likelihood = []

        # estimated allele frequency by EM and LRT
        self.af_by_lrt = {}

        self.depth = {b: 0 for b in self.cmm.BASE}

        quals = np.array(quals)
        self.qual_pvalue = 1.0 - np.exp(self.cmm.MLN10TO10 * quals)
        for i, b in enumerate(bases):
            # Individual likelihood for [A, C, G, T], one sample per row
            if b != 'N' and b[0] not in ['-', '+']:  # ignore all the 'N' bases and indels
                self.ind_allele_likelihood.append([self.qual_pvalue[i]
                                                   if b == t else (1.0 - self.qual_pvalue[i]) / 3
                                                   for t in self.cmm.BASE])
                # record depth for [ACGT]
                if b in self.depth:  # ignore '*'
                    self.depth[b] += 1

        self.ind_allele_likelihood = np.array(self.ind_allele_likelihood)
        self.total_depth = float(sum(self.depth.values()))

        return

    def _set_allele_frequence(self, bases):
        """
        init the base likelihood by bases

        ``bases``: a list like
        """
        total_depth = float(sum([self.depth[b] for b in bases]))
        allele_frequence = np.zeros(len(self.cmm.BASE))  # [A, C, G, T] set to 0.0

        if total_depth > 0:
            for b in bases:
                allele_frequence[self.cmm.BASE2IDX[b]] = self.depth[b] / total_depth

        return np.array(allele_frequence)

    def _f(self, bases, n):
        """
        Calculate population likelihood for all the combination of bases

        Parameters
        ----------

        ``bases``: 1d array like
            A list of bases from [A, C, G, T]

        ``n``: Integer
            The combination number. n must less or equal
            to the length of ``bases``

        Return
        ------

        ``bc``: array=like, combination bases
        ``lr``: Likelihood of ``bc``

        Example
        -------

        >>> import itertools
        >>> bases = ['A', 'C', 'G', 'T']
        >>> bc=[i for i in itertools.combinations(bases,3)]
        >>> bc
        ... [('A', 'C', 'G'), ('A', 'C', 'T'), ('A', 'G', 'T'), ('C', 'G', 'T')]

        """

        bc, lr, bp = [], [], []
        for b in [i for i in itertools.combinations(bases, n)]:
            init_allele_frequecies = self._set_allele_frequence(b)
            if sum(init_allele_frequecies) == 0:
                ## The coverage is empty
                continue

            _, marginal_likelihood, expect_allele_freq = EM(
                np.tile(init_allele_frequecies, (self.ind_allele_likelihood.shape[0], 1)),
                self.ind_allele_likelihood
            )

            bc.append(b)
            lr.append(np.log(marginal_likelihood).sum())  # sum the marginal likelihood
            bp.append(expect_allele_freq)

        return bc, lr, bp

    def lrt(self, specific_base_comb=None):
        """The main function. likelihood ratio test.

        Parameter:
            ``specific_base_comb``: list like
                just calculate the LRT from these specific base combination
        """
        if self.total_depth == 0:
            return

        if specific_base_comb:
            # get effective bases which count frequence >= self.cmm.MINAF
            bases = [b for b in specific_base_comb
                     if self.depth[b] / self.total_depth >= self.cmm.MINAF]
        else:
            # get effective bases which count frequence >= self.cmm.MINAF
            bases = [b for b in self.cmm.BASE
                     if self.depth[b] / self.total_depth >= self.cmm.MINAF]

        if len(bases) == 0:
            return

        # init. Base combination will just be the ``bases`` if specific_base_comb
        bc, lr_null, bp = self._f(bases, len(bases))

        chi_sqrt_value = 0
        base_frq = bp[0]
        lr_alt = lr_null[0]
        for n in range(1, len(bases))[::-1]:

            bc, lr_null, bp = self._f(bases, n)
            lrt_chivalue = 2.0 * (lr_alt - lr_null)
            i_min = np.argmin(lrt_chivalue)

            lr_alt = lr_null[i_min]
            chi_sqrt_value = lrt_chivalue[i_min]
            if chi_sqrt_value < self.cmm.LRT_THRESHOLD:
                # Take the null hypothesis and continue
                bases = bc[i_min]
                base_frq = bp[i_min]
            else:
                # Take the alternate hypothesis
                break

        self._alt_bases = [b for b in bases if b != self._ref_base]
        self.af_by_lrt = {b: '%f' % round(base_frq[self.cmm.BASE2IDX[b]], 6)
                          for b in bases if b != self._ref_base}

        # Todo: improve the calculation method for var_qual
        if len(self._alt_bases):

            r = self.depth[bases[0]] / self.total_depth
            if len(bases) == 1 and self.total_depth > 10 and r > 0.5:
                # mono-allelelic
                self._var_qual = 5000.0

            else:
                chi_prob = chi2.sf(chi_sqrt_value, 1)
                self._var_qual = round(-10 * np.log10(chi_prob), 2) \
                    if chi_prob else 10000.0

            if self._var_qual == 0:
                # _var_qual will been setted as -0.0 instand of 0.0 if it's 0,
                # and I don't know why it so weird!
                self._var_qual = 0.0

        return

    def ref_base(self):
        return self._ref_base

    def alt_bases(self):
        return self._alt_bases

    def var_qual(self):
        return self._var_qual

    def debug(self):
        print(self.ref_base(), self.alt_bases(),
              self.var_qual(), self.depth, self.af_by_lrt)
