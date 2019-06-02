"""Setup file and install script for BaseVar.

Version 1.0.0 (Dec 16, 2018)
Copyright (C) 2018 Shujia Huang <huangshujia9@gmail.com>
"""
import os
from Cython.Distutils import build_ext
from Cython.Build import cythonize

try:
    from setuptools import setup, find_packages, Extension
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages
    from distutils.extension import Extension


DESCRIPTION = "BaseVar: A python software for calling variants from ultra low pass WGS data."
DISTNAME = 'basevar'
MAINTAINER = 'Shujia Huang & Siyang Liu'
MAINTAINER_EMAIL = 'huangshujia9@gmail.com'
URL = 'https://github.com/ShujiaHuang/basevar'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ShujiaHuang/basevar'
VERSION = "0.0.1.3"

#CFLAGS="/home/huangshujia/biosoft/local/include/htslib"
#LDFLAGS="/home/huangshujia/biosoft/local/lib"
BC_INCLUDE_DIR = os.path.split(os.path.realpath(__file__))[0] + "/basevar/caller"
# htslib_dir = os.path.split(os.path.realpath(__file__))[0] + "/basevar/caller/io/htslib"
# cFlags = ["-msse2", "-msse3", "-funroll-loops", "-D_LARGEFILE64_SOURCE", "-D_FILE_OFFSET_BITS=64", "-fPIC"]

CALLER_PRE = 'basevar'
MOD_NAMES = [
    CALLER_PRE + '.io.fasta',
    CALLER_PRE + '.io.bam',
    CALLER_PRE + '.caller.algorithm',
    CALLER_PRE + '.caller.basetype',
]


def make_extension(modname):
    the_cython_file = modname.replace('.', os.path.sep) + '.pyx'
    return Extension(name=modname, sources=[the_cython_file], language='c', include_dirs=[BC_INCLUDE_DIR])


if __name__ == "__main__":
    long_description = os.path.split(os.path.realpath(__file__))[0] + "/README.rst"
    extensions = [make_extension(name) for name in MOD_NAMES]
    
    htslibWrapper = CALLER_PRE+'.io.htslibWrapper'
    the_cython_file = htslibWrapper.replace('.', os.path.sep) + '.pyx'
    extensions.append(Extension(name=htslibWrapper, sources=[the_cython_file], language='c', libraries=['hts']))
    #    Extension(name='htslibWrapper', sources=[(CALLER_PRE+'.io.htslibWrapper').replace(".", os.path.sep)+".pyx"],
    #              language='c', libraries=['hts'], extra_compile_args=cFlags, include_dirs=[htslib_dir]))

    setup(
        name=DISTNAME,
        version=VERSION,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=(open(long_description).read()),
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        packages=find_packages(),
        include_package_data=True,
        ext_modules=cythonize(extensions),
        cmdclass={'build_ext': build_ext},
        install_requires=[
            'Cython==0.29.6',
            'Logbook==1.4.3',
            'numpy==1.15.4',
            'pysam==0.12.0.1',
            'scikit-learn==0.20.2',
            'scipy==1.1.0'
        ],

        # scripts=[],
        entry_points={

            'console_scripts': [
                'basevar = basevar.BaseVar:main'
            ]
        },
        classifiers=[
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.7',
            'License :: OSI Approved :: BSD License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS']
    )
