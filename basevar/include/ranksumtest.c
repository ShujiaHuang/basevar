/*
*   lizilong@bgi.com 201903
*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ranksumtest.h"

static int compare_floats(const void* a, const void* b) {
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
 
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

static double rankR1(double *x, int n1, double *y, int n2) {
    // merge and store the index of sample1
    int ia = n1 - 1;
    int ib = n2 - 1;
    int i, j;
    int *index = malloc(n1 * sizeof(int));
    
    for (i = n1 + n2 - 1; i>=0; i--) {
        if (ia >= 0 && ib < 0) {
            while(ia >= 0){
                index[ia] = ia;
                ia--;
            }
            break;
        }
        if (ia < 0 && ib >= 0) {
            x[i] = y[ib--];
        }
        if (ia >= 0 && ib >= 0) {
            if (x[ia] > y[ib]) {
                index[ia] = i;
                x[i] = x[ia--];
            }else{
                x[i] = y[ib--];
            }
        }
    }
    // average method
    int k = 0, n = 0;
    for (i = 0; i < n1 + n2; i++){
        if (x[i] == x[i+1]) {
            k = k + i + 1;
            n++;
        }else{
            if (k > 0) {
                k = k + i + 1;
                double avg = (double)k / (n + 1);
                for (j = i; j >= i - n; j--) {
                    x[j] = avg;
                }
                k = 0;
                n = 0;
            }else{
                x[i] = i + 1;
                continue;
            }
        }
    }

    //return R1 - the sum of the ranks in Sample1
    double r1 = 0;
    for (i = 0; i < n1; i++) {
        r1 += x[index[i]];
    }

    free(index);
    return r1;
}

double RankSumTest(double *x, int n1, double *y, int n2) {

    int n = n1 + n2;
    double *xx = malloc(n * sizeof(double));
    qsort(x, n1, sizeof(double), compare_floats);
    qsort(y, n2, sizeof(double), compare_floats);
    memcpy(xx, x, n1 * sizeof(double));
    double r1 = rankR1(xx, n1, y, n2);
    double expected = (double)n1 * (n1 + n2 + 1) / 2.0; // fix an unexpected issue 
    double z = (r1 - expected) / sqrt((double)n1*n2*(n1+n2+1)/12.0);

    free(xx);
    return z;
}

