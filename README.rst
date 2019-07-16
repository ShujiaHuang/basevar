BaseVar
=======

Call variants from ultra low-pass WGS data.

Installation
------------

.. code:: bash


    pip install basevar

Prerequisites
-------------

BaseVar requires HTSlib 1.3 or greater. HTSlib can be downloaded from the
`htslib web site <http://www.htslib.org/download/>`_.

To build and install HTSlib, cd into HTSlib source and type `make install`.
This will install HTSlib under `/usr/local/` (see note below). To install HTSlib
in any other directory use `make install prefix=/path/to/dir`.

::

    NOTE: HTSlib should be installed in a standard location (e.g. /usr/local/).

If not installed in a standard location, you will need to set your library paths:

For GNU/Linux
~~~~~~~~~~~~~

.. code:: bash

    export C_INCLUDE_PATH=/path/to/dir/include
    export LIBRARY_PATH=/path/to/dir/lib
    export LD_LIBRARY_PATH=/path/to/dir/lib

Note the `/include` and `/lib` sub-directories. e.g. if you installed HTSlib under `/Users/me/htslib` then set

.. code:: bash

    export C_INCLUDE_PATH=/Users/me/htslib/include
    export LIBRARY_PATH=/Users/me/htslib/lib
    export LD_LIBRARY_PATH=/Users/me/htslib/lib

HTSlib will automatically make the `include` and `lib` directories on install.

For OSX
~~~~~~~

.. code:: bash

    export C_INCLUDE_PATH=/path/to/dir/include
    export LIBRARY_PATH=/path/to/dir/lib
    export DYLD_FALLBACK_LIBRARY_PATH=/path/to/dir/lib


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

