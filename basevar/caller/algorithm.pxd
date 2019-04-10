cdef extern from "include/em.c":
    pass

cdef extern from "include/em.h":
    void em(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood,
            double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon)

cdef extern from "include/ranksumtest.c":
    pass

cdef extern from "include/ranksumtest.h":
    double RankSumTest(double *x, int n1, double *y, int n2)

cdef extern from "include/kfunc.c":
    pass

cdef extern from "include/kfunc.h":
    double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two)
