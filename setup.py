"""Setup file and install script for BaseVar.

Last modify (Jul 16, 2019)
Copyright (C) 2019 Shujia Huang <huangshujia9@gmail.com>
"""
import os
import re

from setuptools import setup, find_packages, Extension
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

from Cython.Distutils import build_ext
from Cython.Build import cythonize

DESCRIPTION = "BaseVar: A python software for calling variants from ultra low pass WGS data."
DISTNAME = 'basevar'
MAINTAINER = 'Shujia Huang & Siyang Liu'
MAINTAINER_EMAIL = 'huangshujia9@gmail.com'
URL = 'https://pypi.org/project/basevar'
DOWNLOAD_URL = 'https://pypi.org/project/basevar'
LICENSE = 'BSD (3-clause)'

ROOT_DIR = os.path.split(os.path.realpath(__file__))[0]


def get_version():
    try:
        f = open(ROOT_DIR + "/basevar/_version.py")
    except EnvironmentError:
        return None

    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            version = mo.group(1)
            return version

    return None


class sdist(_sdist):

    def run(self):
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)


class install(_install):

    def run(self):
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return


BC_INCLUDE_DIR = ROOT_DIR + "/basevar/caller"
TB_INCLUDE_DIR = ROOT_DIR + "/basevar/io/BGZF"

CALLER_PRE = 'basevar'
MOD_NAMES = [
    CALLER_PRE + '.utils',
    CALLER_PRE + '.io.openfile',
    CALLER_PRE + '.io.fasta',
    CALLER_PRE + '.io.bam',
    CALLER_PRE + '.io.read',
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


if __name__ == "__main__":

    # extension for htslib!
    htslibWrapper = CALLER_PRE+'.io.htslibWrapper'
    the_cython_file = htslibWrapper.replace('.', os.path.sep) + '.pyx'
    extensions = [Extension(name=htslibWrapper, sources=[the_cython_file], language='c', libraries=['hts'])]

    algorithm = CALLER_PRE + '.caller.algorithm'
    extensions.append(Extension(name=algorithm,
                                sources=[algorithm.replace('.', os.path.sep) + '.pyx'],
                                language='c',
                                include_dirs=[BC_INCLUDE_DIR],
                                libraries=['hts']))

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
        version=get_version(),
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=(open(ROOT_DIR + "/README.rst").read()),
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        packages=find_packages(),
        include_package_data=True,
        ext_modules=cythonize(extensions),
        cmdclass={'build_ext': build_ext, 'sdist': sdist, 'install': install},
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
                'basevar = basevar.runner:main'
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
