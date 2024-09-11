# BaseVar: Call variants from ultra low-pass WGS data

*BaseVar* was specifically designed to process variant calling from ultra low-depth (<1x) sequencing data, especilly for non-invasive prenatal test (**NIPT**) sequencing data in human genetic studies.

BaseVar has been fully implemented by C++ now. Great improvements were made in the C++ implemetation compare to the [original Python version](https://github.com/ShujiaHuang/basevar/tree/python-version-0.6.1.1). Now, the computing speed of BaseVar is more than 20 times faster than the Python version, and requires much less memory. Generally, each thread (-t/--thread) requires only 3GB to 4GB if -B (--batch-count) option is set to 200, while the Python version need more than 20GB.


## Installation

*BaseVar requires C++17 or higher.*

Build the source codes step-by-step.


### How to install htslib

**1. Download BaseVar from github**

BaseVar is hosted on Github and can be downloaded with the following command:

```bash
$ git clone --recursive https://github.com/ShujiaHuang/basevar.git
```

> WARNING: Please try several times if fail to clone the data causing by 
> the network problem.


**2. Navigate into htslib/htscodecs folder and run the following commands**

After cloing, navigate into the `basevar` folder (cd basevar) and execute the following:

```bash

$ cd htslib/htscodecs
$ autoreconf -i
$ ./configure
$ make

```

**3. Go back to the upper folder and install main htslib by running following commands**

```bash

$ cd htslib
$ autoreconf -i
$ ./configure
$ make

```

**Note**: If you hit something error information looks like the following when you're compilling the `htslib` above, 
you can ignore it and the codes should still work well.

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

Navegate into `bin/` folder (basevar/bin) first and do the following commands:

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

**BaseVar** is under active development. Obtain the newest version by pulling the newest version and compilling again.


## Quick start

### Call variants from several bamfiles

```bash
basevar basetype -R reference.fasta \
    -B 200 -t 4 \
    -I 00alzqq6jw.bam \
    -I 09t3r9n2rg.bam \
    -I 0fkpl1p55b.bam \
    -I 13dg1gvsfk.bam \
    -I 17phildszl.bam \
    -I 1dbpgqt0dq.bam \
    -I 1kyws27hoc.bam \
    --pop-group=sample_group.info \
    --regions=chr11:5246595-5248428,chr17:41197764-41276135 \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz
```

The format of `sample_group.info` could be found [here](tests/data/140k_thalassemia_brca_bam/sample_group.info).


### Or call variants from bamlist

```bash

basevar basetype -R reference.fasta -B 200 -t 4 \
    -L bamfile.list \ 
    --regions=chr11:5246595-5248428,chr17:41197764-41276135 \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz
```

## Citation

Please cite the follow papers if you use BaseVar in your publish projects or papers. 

- Liu et al., [Genomic analyses from non-invasive prenatal testing reveal genetic associations, patterns of viral infections, and chinese population history](https://doi.org/10.1016/j.cell.2018.08.016). *Cell* 175.2 (2018), 347-359

- Liu et al., [Utilizing Non-Invasive Prenatal Test Sequencing Data for Human Genetic Investigation](https://www.biorxiv.org/content/10.1101/2023.12.11.570976v1). BioRxiv (2023)


