BaseVar
=======

Call variants from ultra low-pass WGS data.

Installation
------------

.. code:: bash


    pip install basevar

Quick start
-----------

Call variants from several bamfiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

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

Or call variants from bamlist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash


    basevar basetype -R reference.fasta \
        --regions chr11:5246595-5248428,chr17:41197764-41276135 \
        --batch-count 50 \
        -L bamfile.list \ 
        --output-vcf test.vcf.gz \
        --output-cvg test.cvg.tsv.gz \
        --nCPU 4 && echo "** 5 done **"

