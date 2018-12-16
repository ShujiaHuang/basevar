"""Setup file and install script for BaseVar.

Version 1.0.0 (Dec 16, 2018)
Copyright (C) 2018 Shujia Huang <huangshujia9@gmail.com>
"""
import os

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup


DESCRIPTION = "BaseVar: A python software for calling variants without calling genotype."
DISTNAME = 'basevar'
MAINTAINER = 'Shujia Huang (at) BGI'
MAINTAINER_EMAIL = 'huangshujia9@gmail.com'
URL = 'https://git.bgionline.cn/huangshujia/BaseVar'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://git.bgionline.cn/huangshujia/BaseVar'
VERSION = "1.0.0"


if __name__ == "__main__":

    #long_description = os.path.split(os.path.realpath(__file__))[0] + "/README.rst"
    long_description = os.path.split(os.path.realpath(__file__))[0] + "/README.md"
    requirements_file = os.path.split(os.path.realpath(__file__))[0] + "/requirements.txt"

    requirements = []
    with open(requirements_file) as I:
        for line in I:
            requirements.append(line.strip())

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
        install_requires=requirements,

        scripts=[],
        entry_points={

            'console_scripts': [
                'basevar = basevar.BaseVar:main'
            ]
        },
        classifiers=[
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'License :: OSI Approved :: BSD License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS']
    )
