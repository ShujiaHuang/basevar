# BaseVar: Call variants from ultra low-pass WGS data

*BaseVar* is a specialized tool tailored for variant calling using ultra low-depth (<1x) sequencing data, particularly catering to non-invasive prenatal test (NIPT) data in human genetic studies. Leveraging maximum likelihood and likelihood ratio models, BaseVar accurately identifies polymorphisms at genomic positions and calculates allele frequencies. For in-depth mathematical explanations, refer to the comprehensive documentation available [here](https://doi.org/10.1016/j.cell.2018.08.016).

Now, BaseVar has been fully implemented by C++. BaseVar showcases significant enhancements over its [original Python counterpart](https://github.com/ShujiaHuang/basevar/tree/python-version-0.6.1.1). The C++ implementation delivers a computing speed exceeding 20 times that of the Python version, all while demanding substantially less memory. Typically, each thread (-t/--thread) consumes merely 3GB to 4GB of memory when the -B (--batch-count) option is configured to 200, a stark contrast to the Python version's requirement of over 20GB.


## Citation

Please cite the following paper if you use BaseVar in your published projects or papers. 

> - [1] Liu S, Liu Y, Gu Y, Lin X, Zhu H, Liu H, Xu Z, Cheng S, Lan X, Li L, Huang M, Li H, Nielsen R, Davies RW, Albrechtsen A, Chen GB, Qiu X, Jin X, Huang S. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669 [doi:10.1016/j.xgen.2024.100669](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00288-X)
> - [2] Liu, S., Huang, S., Chen, F., Zhao, L., Yuan, Y., Francis, S.S., Fang, L., Li, Z., Lin, L., Liu, R., et al. (2018). Genomic Analyses from Non-invasive Prenatal Testing Reveal Genetic Associations, Patterns of Viral Infections, and Chinese Population History. *Cell* 175, 347–359.e14. [doi:10.1016/j.cell.2018.08.016](https://doi.org/10.1016/j.cell.2018.08.016).


## Installation

*BaseVar requires C++17 or higher.* Compile `base` from source codes step-by-step.

You can install `basevar` using either of the following two methods.

### Method 1. Install `basevar` by using cmake (Recommend)

This is the simplest way of installing basevar by *cmake*

```bash
$ git clone --recursive https://github.com/ShujiaHuang/basevar.git
$ cd basevar
$ mkdir build
$ cmake ..
$ make 

```

> If you have problems downloading, please try several times.

If everything is smooth, you'll find an exectutable file named `basevar` in `basevar/bin/` folder.

**CAUTION:** If you encounter an error message similar to the following during the compilation of htslib:

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

you can safely disregard it as the code should continue to function properly. To continue the installation processes of basevar, 
you can comment out relevant codes of compiling `htslib` in `CMakeList.txt` by adding "#" at the beginning of each line:

```
add_custom_target(
    Buildhts ALL
    COMMAND autoreconf -i && ./configure && make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
)
```

commnet out the lines:

```
# add_custom_target(
#     Buildhts ALL
#     COMMAND autoreconf -i && ./configure && make
#     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
# )
```

then type `cmake .. && make` again in `build` folder to contiune the following processes.


### Method 2. Manually install processes (Optional)

**1. Download BaseVar from github**

BaseVar is hosted on Github and can be downloaded with the following command:

```bash
$ git clone --recursive https://github.com/ShujiaHuang/basevar.git
```

> **WARNING**: Please try several times if fail to clone the data causing by the network problem.


**2. Navigate into htslib/htscodecs folder and run the following commands**

After cloing, navigate into the `basevar` folder (`cd basevar`) and execute the following:


```bash

$ cd htslib
$ autoreconf -i
$ ./configure
$ make

```

**Note**: If you encounter an error message similar to the following during the compilation of htslib, you can safely disregard it as the code should continue to function properly:

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

Feel free to proceed with your installation tasks despite encountering this error during the compilation process.

**3. Go back to the upper directory and install `basevar` by running the commands below**

Navigate into `bin/` folder (`cd basevar/bin`) first and execute the following commands:

**Manually install in Linux**

```bash
$ cd bin/
$ g++ -O3 -fPIC ../src/*.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto -o basevar


```

**Manually install in MacOS**

```bash
$ cd bin/
$ g++ -O3 -fPIC ../src/*.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o basevar


```

To review each of the parameters, you can type `basevar basetype -h` in the Linux/MacOS Terminal. 

```bash
$ /path/to/basevar basetype -h

About: Call variants and estimate allele frequency by BaseVar.
Usage: basevar basetype [options] <-R Fasta> <--output-vcf> <--output-cvg> [-I input] ...

optional arguments:
  -I, --input=FILE             BAM/CRAM file containing reads.
  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.
  -R, --reference FILE         Input reference fasta file.

  -m, --min-af=float           Setting prior precision of MAF and skip ineffective caller positions,
                               a typical approach involves setting it to min(0.001, 100/x), where x
                               represents the number of input BAM files [min(0.001, 100/x)]. In most
                               cases, users need not be overly concerned about this parameter, as it
                               is generally handled automatically by the program.
  -q, --mapq=INT               Only include reads with mapping quality >= INT. [10]
  -B, --batch-count=INT        INT simples per batchfile. [200]
  -t, --thread=INT             Number of threads. [4]

  -G, --pop-group=FILE         Calculating the allele frequency for specific population.
  -r, --regions=chr:start-end  Skip positions which not in these regions. This parameter could be a list
                               of comma deleimited genome regions(e.g.: chr:start-end) or a file contain
                               the list of regions.
  --output-vcf FILE            Output VCF file.
  --output-cvg FILE            Output position coverage file.

  --filename-has-samplename    If the name of bamfile is something like 'SampleID.xxxx.bam', set this
                               argrument could save a lot of time during get the sample id from BAMfile.
  --smart-rerun                Rerun process by checking batchfiles.
  -h, --help                   Show this help message and exit.
```

This command will provide detailed information about parameters of `basevar`.


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

For stramlinened variant calling across the entire genome, you can use the pipeline generator [**create_pipeline.py**](https://github.com/ShujiaHuang/basevar/blob/master/scripts/create_pipeline.py), which distributes the computational tasks based on the --delta parameter across a specific chromosome defined by the -c parameter.

```bash
$ python create_pipeline.py -R $ref --ref_fai $ref_fai -c chr20 --delta 5000000 -t 20 -L $bamlist -o $outdir > basevar.chr20.sh
```

**BaseVar** is under active development. Obtain the newest version by pulling the newest version and compilling again.

