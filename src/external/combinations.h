// combinations.h, from https://stackoverflow.com/questions/9430568/generating-combinations-in-c

#ifndef __INCLUDE_COMBINATIONS_H__
#define __INCLUDE_COMBINATIONS_H__

#include <vector>

/**
 * @brief This is a recursive method, which you can use on any type. 
 * You can iterate on an instance of Combinations class (e.g. or get() vector with all combinations, 
 * each combination is a vector of objects. This is written in C++11.
 * 
 * modify by Shujia Huang.
 * 
 * @tparam T 
 * 
 */

template<typename T> 
class Combinations {
// Combinations(std::vector<T> s, int m) iterate all Combinations without repetition
// from set s of size m s = {0,1,2,3,4,5} all permuations are: {0, 1, 2}, {0, 1,3}, 
// {0, 1, 4}, {0, 1, 5}, {0, 2, 3}, {0, 2, 4}, {0, 2, 5}, {0, 3, 4}, {0, 3, 5},
// {0, 4, 5}, {1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 3, 4}, {1, 3, 5}, {1, 4, 5}, 
// {2, 3, 4}, {2, 3, 5}, {2, 4, 5}, {3, 4, 5}

public:
    Combinations(std::vector<T> s, int m) : M(m), set(s), partial(std::vector<T>(M)), count(0) {
        N = s.size(); // unsigned long can't be casted to int in initialization

        out = std::vector<std::vector<T>>(comb(N,M), std::vector<T>(M)); // allocate space

        generate(0, N-1, M-1);
    };

    typedef typename std::vector<std::vector<T>>::const_iterator const_iterator;
    typedef typename std::vector<std::vector<T>>::iterator iterator;
    iterator begin() { return out.begin(); }
    iterator end() { return out.end(); }    
    std::vector<std::vector<T>> get() { return out; }

private:
    void generate(int i, int j, int m);
    unsigned long long comb(unsigned long long n, unsigned long long k); // C(n, k) = n! / (n-k)!

    int N, M;
    std::vector<T> set;
    std::vector<T> partial;
    std::vector<std::vector<T>> out;   

    int count;
};

template<typename T> 
void Combinations<T>::generate(int i, int j, int m) {  
    // combination of size m (number of slots) out of set[i..j]
    if (m > 0) { 
        for (int z=i; z<j-m+1; z++) { 
            partial[M-m-1]=set[z]; // add element to permutation
            generate(z+1, j, m-1);
        }
    } else {
        // last position
        for (int z=i; z<j-m+1; z++) { 
            partial[M-m-1] = set[z];
            out[count++] = std::vector<T>(partial); // add to output vector
        }
    }
}

template<typename T> 
unsigned long long Combinations<T>::comb(unsigned long long n, unsigned long long k) {
    // this is from Knuth vol 3

    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

#endif
