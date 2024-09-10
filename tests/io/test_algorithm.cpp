
#include <iostream>
#include <cmath>
#include <string>

#include "algorithm.h"


int main(int argc, char *argv[]) {

    double chi_sqrt_value = -0.1;
    double df = 1;
    double p_value = chi2_test(chi_sqrt_value, df);
    std::cout << "chi2_test(" << chi_sqrt_value << "," << df << ") = " 
              << p_value << " - " << std::isnan(p_value) << "\n";

    // Rank-sum-test
    std::vector<double> sample1 = {1, 2, 3, 10, 7, 3, 3, 4, 5};
    std::vector<double> sample2 = {6, 7, 2, 2, 8, 9, 10};
    double ranksum_test_p = wilcoxon_ranksum_test(sample1, sample2);
    std::cout << "wilcoxon_ranksum_test: " << ranksum_test_p << std::endl;

    // Fisher exact test
    double p = fisher_exact_test(345, 455, 260, 345);
    std::cout << "Fisher exact test: " << p << "\n";
    std::cout << "Fisher exact test: " << fisher_exact_test(8, 4, 4, 9) << "\n";
    std::cout << "Fisher exact test: " << fisher_exact_test(10, 5, 4, 9) << "\n";
    std::cout << "Fisher exact test: " << fisher_exact_test(3, 4, 4, 5) << "\n";
    std::cout << "Fisher exact test: " << fisher_exact_test(1, 1, 1, 1) << "\n";

    // how to test EM?

}
