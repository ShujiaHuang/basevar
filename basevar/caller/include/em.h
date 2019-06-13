#ifndef EM_H
#define EM_H

void em(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood,
        double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon);

#endif
