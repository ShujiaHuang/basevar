// test_combinations.cpp
// compile with: gcc -O3 -Wall -I ../../src -std=c++11 -lstdc++ -o test_combinations test_combinations.cpp
#include <iostream>

#include "external/combinations.h"
#include "utils.h"

struct Bla{
    float x, y, z;
};

int main() {

    std::vector<char> b{'A', 'C', 'G', 'T'};
    std::vector<int> s{0,1,2,3,4,5};
    std::vector<Bla> ss{{1, .4, 5.0},{2, .7, 5.0},{3, .1, 2.0},{4, .66, 99.0}};

    Combinations<char> a(b, 3);
    std::vector<std::vector<char>> ba = a.get();
    // iterate over all combinations
    for (auto x : a) { std::cout << ngslib::join(x, ",") << "\n"; }
    std::cout << "\n1 - "  + ngslib::join(ba[0], "-") << "\n";

    std::cout << "\n\n";

    Combinations<int> c(s,3);
    // iterate over all combinations
    for (auto x : c) { for (auto ii : x) std::cout << ii << ", "; std::cout << "\n"; }

    // or get a vector back
    std::vector<std::vector<int>> z = c.get();  

    std::cout << "\n\n";

    Combinations<Bla> cc(ss, 2);
    // combinations of arbitrary objects
    for (auto x : cc) { for (auto b : x) std::cout << "(" << b.x << ", " << b.y << ", " << b.z << "), "; std::cout << "\n"; }    

}
