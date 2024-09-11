# BaseVar: Call variants from ultra low-pass WGS data

BaseVar has been fully implemented by C++. Great improvements were made in the C++ implemetation compare to the [original Python version](https://github.com/ShujiaHuang/basevar/tree/python-version-0.6.1.1). Now, the computing speed of BaseVar is more than 20 times faster than the Python version, and requires much less memory. Generally, each thread (-t/--thread) requires only 3GB to 4GB if -B (--batch-count) option is set to 200, while the Python version need more than 20GB.


## Installation

Build the source codes step-by-step.


### How to install htslib

**1. Download BaseVar from github**

```bash
$ git clone --recursive https://github.com/ShujiaHuang/basevar.git
```

> WARNING: Please try several times if fail to clone the data causing by 
> the network problem.


**2. Shift to htscodecs directory and run the following commands**

```bash

$ cd htslib/htscodecs
$ autoreconf -i
$ ./configure
$ make

```

**3. Go back to the upper directory and install main htslib by running the commands below**

```bash

$ cd htslib
$ autoreconf -i
$ ./configure
$ make

```

**Note**: If you hit something like the following error information when compile the `htslib` above, we can ignore it
and the codes should still works fine.

```bash
test/test_khash.c: In function 'write_stats_str2int':
test/test_khash.c:53:9: warning: implicit declaration of function 'kh_stats' [-Wimplicit-function-declaration]
   53 |     if (kh_stats(str2int, h, &empty, &deleted, &hist_size, &hist) == 0) {
      |         ^~~~~~~~
test/test_khash.c:53:18: error: 'str2int' undeclared (first use in this function)
   53 |     if (kh_stats(str2int, h, &empty, &deleted, &hist_size, &hist) == 0) {
      |                  ^~~~~~~
test/test_khash.c:53:18: note: each undeclared identifier is reported only once for each function it appears in
make: *** [test/test_khash.o] Error 1
```

**4. Go back to the upper directory and install `basevar` by running the commands below**

**For MacOS**

```bash
$ cd bin/
$ g++ -O3 -fPIC ../src/main.cpp ../src/basetype.h ../src/basetype.cpp ../src/basetype_caller.cpp ../src/utils.cpp ../src/fasta.cpp ../src/bam_header.cpp ../src/bam.cpp ../src/bam_record.cpp ../src/basetype_utils.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o basevar

```

**For Linux**

```bash
$ cd bin/
$ g++ -O3 -fPIC ../src/main.cpp ../src/basetype.h ../src/basetype.cpp ../src/basetype_caller.cpp ../src/utils.cpp ../src/fasta.cpp ../src/bam_header.cpp ../src/bam.cpp ../src/bam_record.cpp ../src/basetype_utils.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto -o basevar

```

## Quick start

### Call variants from several bamfiles

```bash
basevar basetype -R reference.fasta \
    -I 00alzqq6jw.bam \
    -I 09t3r9n2rg.bam \
    -I 0fkpl1p55b.bam \
    -I 13dg1gvsfk.bam \
    -I 17phildszl.bam \
    -I 1dbpgqt0dq.bam \
    -I 1kyws27hoc.bam \
    --pop-group=sample_group.info \
    --regions=chr11:5246595-5248428,chr17:41197764-41276135 \
    -B 200 -t 4 \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz
```

The format of `sample_group.info` could be found [here](tests/data/140k_thalassemia_brca_bam/sample_group.info).


### Or call variants from bamlist

```bash

basevar basetype -R reference.fasta -L bamfile.list \ 
    --regions=chr11:5246595-5248428,chr17:41197764-41276135 \
    -B 200 -t 4 \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz
```

