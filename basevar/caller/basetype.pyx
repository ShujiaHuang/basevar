# cython: profile=True
"""
This module contain functions of LRT and Base genotype.
"""
import itertools  # Use the combinations function

from scipy.stats.distributions import chi2

from basevar.utils import CommonParameter
from basevar.caller.algorithm cimport EM


cdef class BaseTuple:
    def __cinit__(self, int combination_num, int base_num, int base_type_num):

        # initial the value
        self.combination_num = combination_num
        self.base_num = base_num
        self.base_comb_tuple = <char**>(calloc(combination_num, sizeof(char*)))
        assert self.base_comb_tuple != NULL, "Could not allocate memory for self.base_comb_tuple in BaseTuple."

        self.alleles_freq_list = <double**>(calloc(combination_num, sizeof(double*)))
        assert self.alleles_freq_list != NULL, (
            "Could not allocate memory for alleles_freq_list in BaseTuple.")

        cdef int i = 0
        cdef int j = 0
        for i in range(combination_num):
            self.base_comb_tuple[i] = <char*>(malloc(base_num * sizeof(char)))
            self.alleles_freq_list[i] = <double*>(calloc(base_type_num, sizeof(double))) # For A/C/G/T

        self.sum_marginal_likelihood = <double*>(calloc(combination_num, sizeof(double)))
        assert self.sum_marginal_likelihood != NULL, (
            "Could not allocate memory for sum_marginal_likelihood in BaseTuple.")

    def __dealloc__(self):
        self.destroy()

    cdef void destroy(self):
        """Free memory"""
        cdef int index
        if self.base_comb_tuple != NULL:
            for index in range(self.combination_num):
                if self.base_comb_tuple[index] != NULL:
                    free(self.base_comb_tuple[index])
                    self.base_comb_tuple[index] = NULL

            free(self.base_comb_tuple)
            self.base_comb_tuple = NULL

        if self.alleles_freq_list != NULL:
            for index in range(self.combination_num):
                if self.alleles_freq_list[index] != NULL:
                    free(self.alleles_freq_list[index])
                    self.alleles_freq_list[index] = NULL

            free(self.alleles_freq_list)
            self.alleles_freq_list = NULL

        if self.sum_marginal_likelihood != NULL:
            free(self.sum_marginal_likelihood)
            self.sum_marginal_likelihood = NULL


cdef class BaseType:

    def __cinit__(self):
        # do nothings
        pass

    cdef void cinit (self, bytes ref_base, char **bases, int *quals, int total_sample_size, float min_af):
        """ Iinitial all the data here.
        
        A class for calculate the base probability

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
        """
        self._ref_base = ref_base  # ref_base must be upper() before pass to this class.
        self._alt_bases = None
        self._var_qual = 0.0  # init the variant quality
        self.min_af = min_af
        self.depth = {b: 0 for b in CommonParameter.BASE}
        self.base_type_num = len(CommonParameter.BASE)

        # how many individual is in good base
        cdef int i = 0
        self.good_individual_num = 0
        for i in range(total_sample_size):
            if bases[i][0] not in ['N', '-', '+']:
                self.good_individual_num += 1

        # qual_pvalue has to be for all, because we'll use this outside.
        self.qual_pvalue = <double*>(calloc(total_sample_size, sizeof(double)))
        assert self.qual_pvalue != NULL, "Could not allocate memory for qual_pvalue in BaseType"

        for i in range(total_sample_size):
            self.qual_pvalue[i] = 1.0 - exp(CommonParameter.MLN10TO10 * quals[i])

        # A big 1-D array
        self.ind_allele_likelihood = <double*>(calloc(self.good_individual_num * self.base_type_num, sizeof(double)))
        assert self.ind_allele_likelihood != NULL, "Could not allocate memory for ind_allele_likelihood in BaseType"

        # set allele likelihood for each individual and get depth
        self._set_init_ind_allele_likelihood(bases, CommonParameter.BASE, total_sample_size)
        self.total_depth = float(sum(self.depth.values()))

        # estimated allele frequency by EM and LRT
        self.af_by_lrt = {}

        return

    def __dealloc__(self):
        """
        Free memory
        """
        if self.ind_allele_likelihood != NULL:
            free(self.ind_allele_likelihood)

        if self.qual_pvalue != NULL:
            free(self.qual_pvalue)

    cdef void _set_init_ind_allele_likelihood(self, char **ind_bases, list base_element, int total_individual_num):

        cdef int i = 0
        cdef int j = 0
        cdef int k = 0
        for i in range(total_individual_num):

            # Individual likelihood for [A, C, G, T], one sample per row
            # ignore all the 'N' bases and indels.
            if ind_bases[i][0] not in ['N', '-', '+']:

                # Just set allele likelihood for good individual
                for k in range(self.base_type_num):
                    if ind_bases[i] == base_element[k]:
                        self.ind_allele_likelihood[j * self.base_type_num + k] = self.qual_pvalue[i]
                    else:
                        self.ind_allele_likelihood[j * self.base_type_num + k] = (1.0 - self.qual_pvalue[i])/3

                # iteration good individual
                j += 1

                # record coverage for [ACGT]
                if ind_bases[i] in self.depth:
                    self.depth[ind_bases[i]] += 1
        return

    cdef double* _set_allele_frequence(self, tuple bases):
        """
        init the base likelihood by bases

        ``bases``: a list like
        """
        # initial [A, C, G, T] to be 0.0
        cdef double* allele_frequence = <double*>(calloc(self.base_type_num, sizeof(double)))
        assert allele_frequence != NULL, (
            "Could not allocate memory for allele_frequence in BaseType._set_allele_frequence")

        cdef bytes b
        if self.total_depth > 0:
            for b in bases:
                allele_frequence[CommonParameter.BASE2IDX[b]] = self.depth[b] / self.total_depth

        return allele_frequence

    cdef BaseTuple _f(self, list bases, int n):
        """
        Calculate population likelihood for all the combination of bases

        Parameters
        ----------
        ``bases``: 1d array like
            A list of bases from [A, C, G, T]

        ``n``: Integer
            The combination number. n must less or equal
            to the length of ``bases``

        Example
        -------

        >>> import itertools
        >>> bases = ['A', 'C', 'G', 'T']
        >>> bc=[i for i in itertools.combinations(bases,3)]
        >>> bc
        ... [('A', 'C', 'G'), ('A', 'C', 'T'), ('A', 'G', 'T'), ('C', 'G', 'T')]

        """
        cdef double* init_allele_frequecies = NULL
        cdef double* marginal_likelihood = NULL
        cdef double* expect_allele_freq = NULL

        cdef list base_combs_tuple = [x for x in itertools.combinations(bases, n)]
        cdef int comb_num = len(base_combs_tuple)

        cdef BaseTuple base_tuple = BaseTuple(comb_num, n, self.base_type_num)
        cdef int bi = 0
        cdef int i = 0
        for i in range(comb_num):

            # initial the allele frequencies of [A, C, G, T]
            init_allele_frequecies = self._set_allele_frequence(base_combs_tuple[i])
            if self.sum_likelihood(init_allele_frequecies, self.base_type_num, False) == 0:
                free(init_allele_frequecies)
                continue

            # reset every time
            marginal_likelihood = <double*>(calloc(self.good_individual_num, sizeof(double)))
            expect_allele_freq = <double*>(calloc(self.base_type_num, sizeof(double)))

            EM(init_allele_frequecies,
               self.ind_allele_likelihood,
               marginal_likelihood, # update every loop
               expect_allele_freq,  # update every loop
               self.good_individual_num,
               self.base_type_num,
               100,  # EM iter_num
               0.001) # EM epsilon

            for bi in range(n):
                # each element is single base
                base_tuple.base_comb_tuple[i][bi] = ord(base_combs_tuple[i][bi])

            # Todo: Should we use log10 function instead of using log or not? check it carefully!
            # sum the marginal likelihood
            base_tuple.sum_marginal_likelihood[i] = self.sum_likelihood(
                marginal_likelihood, self.good_individual_num, True)

            base_tuple.alleles_freq_list[i] = expect_allele_freq

            free(marginal_likelihood)
            free(init_allele_frequecies)

        expect_allele_freq = NULL
        return base_tuple

    cdef double sum_likelihood(self, double* data, int num, bint is_log):
        cdef double s = 0.0
        cdef int i = 0
        for i in range(num):
            if is_log:
                s += log(data[i])
            else:
                s += data[i]

        # a double-type value
        return s

    cdef bint lrt(self, list specific_base_comb):
        """The main function. likelihood ratio test.

        Parameter:
            ``specific_base_comb``: list like
                Calculating LRT for specific base combination
        """
        if self.total_depth == 0:
            return False

        cdef list bases = []
        if specific_base_comb:
            bases = [b for b in specific_base_comb if self.depth[b] / self.total_depth >= self.min_af]
        else:
            bases = [b for b in CommonParameter.BASE if self.depth[b] / self.total_depth >= self.min_af]

        cdef int bases_num = len(bases)
        if bases_num == 0 or (bases_num == 1 and bases[0] == self._ref_base): # no base or it's reference base.
            return False

        # init. Base combination will just be the ``bases`` if specific_base_comb
        cdef BaseTuple the_base_tuple = self._f(bases, bases_num)
        cdef double lr_alt = the_base_tuple.sum_marginal_likelihood[0]
        cdef double* base_frq = <double*>(calloc(self.base_type_num, sizeof(double)))
        memcpy(base_frq, the_base_tuple.alleles_freq_list[0], self.base_type_num * sizeof(double))

        cdef double chi_sqrt_value = 0
        cdef double* lrt_chi_value = NULL
        cdef int n
        cdef int i_min
        for n in range(1, len(bases))[::-1]:  # From complex to simplicity

            the_base_tuple = self._f(bases, n)
            if lrt_chi_value != NULL:
                free(lrt_chi_value)

            lrt_chi_value = self.calculate_chivalue(lr_alt,
                                                    the_base_tuple.sum_marginal_likelihood,
                                                    the_base_tuple.combination_num)

            i_min = self.find_argmin(lrt_chi_value, the_base_tuple.combination_num)
            lr_alt = the_base_tuple.sum_marginal_likelihood[i_min]
            chi_sqrt_value = lrt_chi_value[i_min]

            # Take the null hypothesis and continue
            if chi_sqrt_value < CommonParameter.LRT_THRESHOLD:
                memcpy(base_frq, the_base_tuple.alleles_freq_list[i_min], self.base_type_num * sizeof(double))
                bases = [chr(the_base_tuple.base_comb_tuple[i_min][j]) for j in range(the_base_tuple.base_num)]

            # Take the alternate hypothesis
            else:
                break

        # clear the_base_tuple
        the_base_tuple.destroy()

        if lrt_chi_value != NULL:
            free(lrt_chi_value)

        self._alt_bases = [b for b in bases if b != self._ref_base]
        self.af_by_lrt = {b:"%.6f" % base_frq[CommonParameter.BASE2IDX[b]] for b in self._alt_bases}

        cdef bint is_variant = False
        # Todo: improve the calculation method for var_qual
        cdef double r
        if len(self._alt_bases):

            is_variant = True

            r = self.depth[bases[0]] / self.total_depth
            if len(bases) == 1 and self.total_depth > 10 and r > 0.5:
                # mono-allelelic
                self._var_qual = 5000.0

            else:
                chi_prob = chi2.sf(chi_sqrt_value, 1)
                self._var_qual = round(-10 * log10(chi_prob)) if chi_prob > 0 else 10000.0

            if self._var_qual == 0:
                self._var_qual = 0.0

        return is_variant

    cdef double* calculate_chivalue(self, double lr_alt, double* lr_null, int comb_num):

        cdef double* chi_value = <double*>(calloc(comb_num, sizeof(double)))
        cdef int i = 0
        for i in range(comb_num):
            chi_value[i] = 2 * (lr_alt - lr_null[i])

        return chi_value

    cdef int find_argmin(self, double* data, int comb_num):
        """Return indices of the minimum values along the given axis of `data`. """
        cdef int i = 0
        cdef int index = 0
        cdef double min_value = data[index]
        for i in range(1, comb_num):
            if data[i] < min_value:
                index = i
                min_value = data[i]

        return index

    property alt_bases:
        """Return the list of variants"""
        def __get__(self):
            return self._alt_bases

    property var_qual:
        """Return the quality score for the variant site"""
        def __get__(self):
            # A double value
            return self._var_qual
