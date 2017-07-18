import itertools   # Use the combinations function
import numpy as np

from algorithm import EM, strand_bias

class CommonParameter(object):
    """
    defined some globle common parameters
    """
    def __init__(self):
        self.LRT_THRESHOLD = 24   ## 24 corresponding to a chi-pvalue of 10^-6
        self.QUAL_THRESHOLD = 60  ## -10 * lg(10^-6)
        self.MLN10TO10 = -0.23025850929940458 # -np.log(10)/10
        self.BASE = ['A', 'C', 'G', 'T']
        self.BASE2IDX = {'A':0, 'C':1, 'G':2, 'T':3}
        self.debug = False
        self.verbose = False
        self.MINAF = 0.001  # The effective base freqence threshold for 140k sample size
        self.scan_indel = True

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
        self._var_qual = 0 # init the variant quality
        self._ref_base = ref_base
        self.cmm = cmm

        #The prior probability for echo base
        self.prior_prob = []
        self.eaf = {} ## estimated allele frequency by EM
        self.depth = {b:0 for b in self.cmm.BASE}

        quals = np.array(quals)
        self.qual_pvalue = 1.0 - np.exp(self.cmm.MLN10TO10 * quals)
        for i, b in enumerate(bases):
            ## Individual per row for [A, C, G, T] and give a prior probability
            if b != 'N':  # ignore all the 'N' base sample
                self.prior_prob.append([self.qual_pvalue[i]
                                        if b == t else (1.0-self.qual_pvalue[i])/3
                                        for t in self.cmm.BASE])
                # record depth for [ACGT]
                if b in self.depth: # ignore '*'
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
        print(self.ref_base(), self.alt_bases(),
              self.var_qual(), self.depth, self.eaf)

    def _set_base_likelihood(self, bases):
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
        bc = []
        lr = []
        bp = []

        for b in [i for i in itertools.combinations(bases, n)]:

            init_likelihood = self._set_base_likelihood(b)
            if sum(init_likelihood) == 0: continue  ## No covert

            _, pop_likelihood, base_expected_prob = EM(
                self.prior_prob,
                np.tile(init_likelihood, (self.prior_prob.shape[0],1)))

            bc.append(b)
            lr.append(np.log(pop_likelihood).sum()) # sum the marginal prob
            bp.append(base_expected_prob)

        return bc, lr, bp

    def lrt(self):
        """The main function.
        likelihood ratio test.
        """
        if self.total_depth == 0: return

        # get effective bases which count frequence > self.cmm.MINAF
        bases = [b for b in self.cmm.BASE
                 if self.depth[b]/self.total_depth > self.cmm.MINAF]

        # init
        _, lr_null, bp = self._f(bases, len(bases))

        base_frq = bp[0]
        lr_alt = lr_null[0]

        chi_sqrt_value = 0
        for n in range(1, len(bases))[::-1]:

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
        self.eaf = {b : '%f' % round(base_frq[self.cmm.BASE2IDX[b]], 6)
                    for b in bases if b != self._ref_base}

        # calculate the variant quality
        # Todo: improve the calculation method for var_qual
        if len(self._alt_bases):

            if len(bases) == 1 and self.depth[bases[0]] / self.total_depth > 0.5:
                # mono-allelelic
                self._var_qual = 5000

            else:
                from scipy.stats import chisqprob as sp_chisqprob
                chi_prob = sp_chisqprob(chi_sqrt_value, 1)
                self._var_qual = round(-10 * np.log10(chi_prob), 2) \
                    if chi_prob else 10000

        return
