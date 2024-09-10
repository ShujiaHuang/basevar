# BaseVar
Call variants from ultra low-pass WGS data.

## Installation

Install the released version by `pip`:

```bash
pip install basevar
```

You may instead want to install the development version from github, by running:

```bash
pip install git+git://github.com/ShujiaHuang/basevar.git#egg=basevar
```


## Quick start

### Call variants from several bamfiles

```bash
basevar basetype -R reference.fasta \
    --regions chr11:5246595-5248428,chr17:41197764-41276135 \
    --batch-count 50 \
    -I 00alzqq6jw.bam \
    -I 09t3r9n2rg.bam \
    -I 0fkpl1p55b.bam \
    -I 13dg1gvsfk.bam \
    -I 17phildszl.bam \
    -I 1dbpgqt0dq.bam \
    -I 1kyws27hoc.bam \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz \
    --nCPU 4 && echo "** 5 done **"
```

### Or call variants from bamlist
```bash

basevar basetype -R reference.fasta \
    --regions chr11:5246595-5248428,chr17:41197764-41276135 \
    --batch-count 50 \
    -L bamfile.list \ 
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz \
    --nCPU 4 && echo "** 5 done **"

```