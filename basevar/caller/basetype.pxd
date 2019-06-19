"""Header for basetype.pyx
"""
cdef extern from "stdlib.h":
    void *malloc(size_t)
    void *calloc(size_t, size_t)
    void *memcpy(void *dst, void *src, size_t length)
    void free(void *)

cdef extern from "math.h":
    double exp(double)
    double round(double)
    double log(double)
    double log10(double)


cdef class BaseTuple:
    cdef int combination_num
    cdef int base_num
    cdef char** base_comb_tuple  # => bc
    cdef double* sum_marginal_likelihood  # => lr
    cdef double** alleles_freq_list  # => bp
    cdef void destory_tuple(self)


cdef class BaseType:
    cdef int good_individual_num
    cdef int base_type_num
    cdef bytes _ref_base
    cdef list _alt_bases
    cdef double total_depth
    cdef double _var_qual
    cdef float min_af
    cdef double* ind_allele_likelihood
    cdef double* qual_pvalue
    cdef dict af_by_lrt
    cdef dict depth

    cdef bint lrt(self, list specific_base_comb)
    cdef void _set_init_ind_allele_likelihood(self, list bases, list base_element)
    cdef double* _set_allele_frequence(self, tuple bases)
    cdef double sum_likelihood(self, double* data, int num, bint is_log)
    cdef BaseTuple _f(self, list bases, int n)
    cdef double* calculate_chivalue(self, double lr_alt, double* lr_null, int comb_num)
    cdef int find_argmin(self, double* data, int comb_num)
