#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

#include <array>
#include <vector>
#include "utils.h"

static const double MLN10TO10 = -0.23025850929940458;  // 换底，把 10 换为 e 为底：log(10)/10;

int main(int argc, char *argv[]) {

std::vector<double> tt{1,2,3,4};
std::cout << ngslib::join(tt, ",") <<"\n";
tt = std::vector<double> (3, 0); 
std::cout << ngslib::join(tt, ",") <<"\n";
tt.clear();

std::vector<char> cc = {'z', 'c', 'y', 'a', 'G', 'C', 'T', 'A'};
std::cout << "Before sort: " << ngslib::join(cc, " - ") <<"\n";

std::sort(cc.begin(), cc.end());
std::cout << "After sort : " << ngslib::join(cc, " - ") <<"\n";

    std::array<int, 7> numbers{ 2, 4, 8, 0, 6, -1, 3};
	int minIndex = ngslib::argmin(numbers.begin(), numbers.end());
	std::cout << "MinIndex: " << minIndex << '\n';
	std::vector<float> prices = { 12.5f, 8.9f, 100.0f, 24.5f, 30.0f };
	float maxIndex = ngslib::argmax(prices.begin(), prices.end());
	std::cout << "MaxIndex: " <<  maxIndex << '\n';

exit(1);

    std::cout << "phred 'a': " << exp(('a'-33) * MLN10TO10) << "\n";

    std::string str = "AbcG";
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    std::cout << "Toupper case: 'AbcG' => " << str << "\n";

    std::vector<std::string> a, b;
    a.push_back("Hello");
    a.push_back("World!");
    std::string s = "Hello world. There are two needles in this haystack with needles.";
    std::cout << "Find first 'are': " << s.find("are") << " - " << (s.find("are", 100) ==std::string::npos) << "\n";
    std::cout << s.size() << " - " << s.length() << "\n";

    const char *ts = "##Sample IDs=ERS225193";
    std::vector<std::string> h; ngslib::split(ts, h, "=");
    std::cout << "##Sample IDs=ERS225193 : " << h[0] << " - " << h[1] << "\n" << ngslib::join(h, "/") << "\n";
    
    ngslib::split(". I am hello\tworld,w", b, ".");
    std::cout << "'. I am thello\\tworld,w' : " << b.size() << " : " << ngslib::join(b, "--") << std::endl;

    ngslib::split("hello, world :: ww", b, ":");
    std::cout << "'hello, world :: ww'    : " << b.size() << " " << ngslib::join(b, "--") << std::endl;

    ngslib::split(",hello,world,,ww,", b, ",");
    std::cout << "',hello,world,,ww,'  : " << b.size() << " " << ngslib::join(b, "--") << std::endl;

    ngslib::split("hello::world", b, ">>", false);
    ngslib::split("hello::world", b, "h", false);
    std::cout << b[0] << ": " << ngslib::join(b, "--") << std::endl;

    std::vector<int> c = {1,2,3,4};
    std::cout << ngslib::join(c) << std::endl;
    std::cout << ngslib::join(c, ":") << std::endl;

    std::vector<float> d = {0.1,2.5,3,4.0};
    std::cout << ngslib::join(d, " - ") << std::endl;

    ngslib::split("1,2,3,4", c, ",", false);
    ngslib::split(",1.1,2,3,4", d, ",", true);  // The first element is empty and will be initial by 0.
    std::cout << ngslib::join(c, " ") << std::endl;
    std::cout << ngslib::join(d, " ") << std::endl;

    std::cout << "Get filename from Linux path: " << ngslib::basename("some/path/file.ext") << "\t" 
              << ngslib::remove_filename_extension(ngslib::basename("some/path/file.ext")) << "\n";
    std::cout << "Get filename from window path " << ngslib::basename("C:\\\\MyDirectory\\\\MyFile.bat") << "\t" 
              << ngslib::remove_filename_extension("C:\\\\MyDirectory\\\\MyFile.bat") << "\n";

    std::cout << "Get dirname from Linux path: " << ngslib::dirname("some/path/file.ext") << "\t" 
              << ngslib::dirname("./a") << "\n";
    std::cout << "Get dirname from window path " << ngslib::dirname("C:\\\\MyDirectory\\\\MyFile.bat") << "\n";

    std::cout << "is_readable: " << ngslib::is_readable("some") << " - test_threadpool - " << ngslib::is_readable("test_bamrecord") << "\n";
    std::cout << ngslib::safe_mkdir("a/b/c")  << "\n";
    std::cout << ngslib::safe_remove("a/b/c") << "\n";
    std::cout << ngslib::safe_remove("t_t_")  << "\n";

    std::string lmf = ngslib::get_last_modification_file("../data/");
    std::cout << "-- LMF: " << lmf << "\n";
    std::cout << "-- LMF a/b: " << ngslib::get_last_modification_file("a/b") << "\n";
    std::cout << "-- LMF ../data/bam100: " << ngslib::get_last_modification_file("../data/bam100") << "\n";

    std::vector<int> vc = ngslib::vector_slicing(c, 1, 3);  // [start, end)
    std::cout << "vector slicing `vc` size: " << vc.size() << "\t" << ngslib::join(vc, "-") << "\n";

    std::cout << "Get Linux absolute path: " << ngslib::abspath("some/path/file.ext") << "\n";
    std::cout << "Get linux absolute path: " << ngslib::abspath("/file.ext") << "\n";

    return 0;
}