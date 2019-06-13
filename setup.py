"""Setup file and install script for BaseVar.

Version 1.0.0 (Dec 16, 2018)
Copyright (C) 2018 Shujia Huang <huangshujia9@gmail.com>
"""
import os
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

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

BC_INCLUDE_DIR = os.path.split(os.path.realpath(__file__))[0] + "/basevar/caller"
TB_INCLUDE_DIR = os.path.split(os.path.realpath(__file__))[0] + "/basevar/io/BGZF"

CALLER_PRE = 'basevar'
MOD_NAMES = [
    CALLER_PRE + '.utils',
    CALLER_PRE + '.io.openfile',
    CALLER_PRE + '.io.fasta',
    CALLER_PRE + '.io.bam',
    CALLER_PRE + '.io.read',
    # CALLER_PRE + '.caller.algorithm',
    CALLER_PRE + '.caller.basetype',
    CALLER_PRE + '.caller.batch',
    CALLER_PRE + '.caller.batchcaller',
    CALLER_PRE + '.caller.basetypebam',
    CALLER_PRE + '.caller.basetypeprocess',
    CALLER_PRE + '.caller.executor',
]


def make_extension(modname):
    the_cython_file = modname.replace('.', os.path.sep) + '.pyx'
    return Extension(name=modname, sources=[the_cython_file], language='c')
    # return Extension(name=modname, sources=[the_cython_file], language='c', include_dirs=[BC_INCLUDE_DIR])


if __name__ == "__main__":
    long_description = os.path.split(os.path.realpath(__file__))[0] + "/README.rst"

    # extension for htslib!
    htslibWrapper = CALLER_PRE+'.io.htslibWrapper'
    the_cython_file = htslibWrapper.replace('.', os.path.sep) + '.pyx'
    extensions = [Extension(name=htslibWrapper, sources=[the_cython_file], language='c', libraries=['hts'])]

    algorithm = CALLER_PRE + '.caller.algorithm'
    extensions.append(Extension(name=algorithm,
                                sources=[algorithm.replace('.', os.path.sep) + '.pyx'],
                                language='c',
                                include_dirs=[BC_INCLUDE_DIR]))

    # extension for bgzip and tabix
    # the_tabix_pre = CALLER_PRE+'.io.BGZF.tabix'
    # the_tabix_dir = the_tabix_pre.replace('.', os.path.sep)
    # the_pysam_pre = CALLER_PRE+'.io.BGZF.pysam'
    # the_pysam_dir = the_pysam_pre.replace('.', os.path.sep)
    # tabix_flags = ["-Wno-incompatible-pointer-types-discards-qualifiers","-Wno-unused-function","-Wno-unneeded-internal-declaration"]
    # tabproxies_flags = ["-Wno-unused-function"]
    # extensions.append(Extension(name=the_pysam_pre+".ctabix",
    #                             sources=[the_pysam_dir+"/ctabix.pyx"] + [the_pysam_dir+"/tabix_util.c"] + glob.glob(the_tabix_dir+"/*.pysam.c"),
    #                             include_dirs=[TB_INCLUDE_DIR+"/tabix", TB_INCLUDE_DIR+"/pysam"],
    #                             libraries=["z"],
    #                             language="c",
    #                             extra_compile_args=tabix_flags
    #                             ))
    #
    # extensions.append(Extension(name=the_pysam_pre+".TabProxies",
    #                             sources=[the_pysam_dir+"/TabProxies.pyx"],
    #                             libraries=["z"],
    #                             language="c",
    #                             extra_compile_args=tabproxies_flags
    #                             ))

    # other extensions
    extensions += [make_extension(name) for name in MOD_NAMES]

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
