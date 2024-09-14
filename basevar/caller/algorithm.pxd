
cdef extern from "stdlib.h":
    void *malloc(size_t)
    void *memcpy(void *dst, void *src, size_t length)
    void free(void *)

cdef extern from "c/em.c":
    pass

cdef extern from "c/em.h":
    void em(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood,
            double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon)

cdef extern from "c/ranksumtest.c":
    pass

cdef extern from "c/ranksumtest.h":
    double RankSumTest(double *x, int n1, double *y, int n2)

cdef void EM(double* init_allele_freq,
             double* ind_allele_likelihood,
             double* marginal_likelihood,
             double* expect_allele_prob,
             int nsample,
             int ntype,
             int iter_num,
             double epsilon)

cdef tuple strand_bias(bytes ref_base, list alt_bases, char **bases, char *strands, int size)
cdef double ref_vs_alt_ranksumtest(bytes ref_base, list alt_base, char **bases, int *info, int data_size)

