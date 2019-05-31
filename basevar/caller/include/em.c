/*
*   lizilong@bgi.com 201903
*/
#include <math.h>
#include <stdlib.h>
#include "em.h"

static void singleEM(double *allele_freq, double *ind_allele_likelihood, double *marginal_likelihood,
                     double *expect_allele_prob, int nsample, int ntype) {
    double *likelihood = (double *) calloc(ntype, sizeof(double));
    double *ind_allele_prob = (double *) calloc(ntype * nsample, sizeof(double)); // covert ind_allele_prob to col major
    // step E 
    int i, j;
    for(i=0; i<nsample; ++i){
        for(j=0; j<ntype; ++j){
            likelihood[j] = allele_freq[j] * ind_allele_likelihood[i * ntype + j];
            marginal_likelihood[i] += likelihood[j];
        }
        for(j=0; j<ntype; ++j){
            /* col major may be fast for step M
             * need to deal with marginal_likelihood[i] is close to zero */
            ind_allele_prob[j * nsample + i] = likelihood[j] / marginal_likelihood[i]; 
        }
    }
    free(likelihood);

    // step M
    for(j=0; j<ntype; ++j){
        for(i=0; i<nsample; ++i){
            expect_allele_prob[j] += ind_allele_prob[j * nsample + i]; 
        }
        expect_allele_prob[j] = expect_allele_prob[j] / nsample;
    }
    free(ind_allele_prob);
    
    return;
}

static void update_allele_freq(double *allele_freq, double *expect_allele_prob, int ntype) {

    int j;
    for(j=0; j<ntype; ++j){
        allele_freq[j] = expect_allele_prob[j];
        expect_allele_prob[j] = 0.0;
    }
    return;
}

static double delta_bylog(double *bf, double *af, int n) {
    double delta = 0.0;
    int i;
    for(i=0; i<n; ++i){
        // need to deal with log(0) == inf;
        delta += fabs(log(af[i]) - log(bf[i]));
        bf[i] = af[i];
        af[i] = 0.0;
    }
    return delta;
}

void em(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood,
        double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon) {

    double *af_marginal_likelihood = (double *) calloc(nsample, sizeof(double));
    double *allele_freq = (double *) malloc(ntype * sizeof(double));
    double delta;
    int i, j;

    /*
     copy allele_freq in case that init_allele_freq be modified; 
    */
    int i, j;
    for(j = 0; j < ntype; ++j){
	    allele_freq[j] = init_allele_freq[j];
    }
    singleEM(allele_freq, ind_allele_likelihood, marginal_likelihood, expect_allele_prob, nsample, ntype);

    for(i=0; i<iter_num; ++i){
        update_allele_freq(allele_freq, expect_allele_prob, ntype);
        singleEM(allele_freq, ind_allele_likelihood, af_marginal_likelihood, expect_allele_prob, nsample, ntype);
        delta = delta_bylog(marginal_likelihood, af_marginal_likelihood, nsample);
        if(delta < epsilon){
            break;
        }
    }

    free(af_marginal_likelihood);
    free(allele_freq);

    return;
}

