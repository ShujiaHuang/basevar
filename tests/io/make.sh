g++ -O3 -fPIC test_fasta.cpp ../../src/fasta.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_fasta && ./test_fasta

g++ -O3 -fPIC test_bamheader.cpp ../../src/bam_header.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bamheader && ./test_bamheader

g++ -O3 -fPIC test_bamrecord.cpp ../../src/bam_record.cpp ../../src/fasta.cpp ../../src/bam_header.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bamrecord && ./test_bamrecord

g++ -O3 -fPIC test_bam.cpp ../../src/bam.cpp ../../src/bam_header.cpp ../../src/bam_record.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bam && ./test_bam

g++ -O3 -fPIC test_cmdline.cpp -I ../../src -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_cmdline && ./test_cmdline

g++ -O3 -fPIC test_threadpool.cpp -I ../../src -o test_threadpool

g++ -O3 -fPIC test_utils.cpp ../../src/utils.cpp -I ../../src -o test_utils
g++ -O3 -fPIC test_algorithm.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -o test_algorithm && ./test_algorithm


gcc -O3 -Wall -I ../../src -std=c++11 -lstdc++ -o test_combinations test_combinations.cpp


g++ -O3 -fPIC test_algorithm.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -o test_algorithm && ./test_algorithm

g++ -O3 -fPIC test_vcfconcat.cpp ../../htslib/libhts.a ../../src/external/vcfconcat.c -I ../../src -I ../../htslib -o test_vcfconcat

