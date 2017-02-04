# Copyright (C) 2016 Shujia Huang <huangshujia9@gmail.com>
import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': ('BaseVar: A python software to call varaint without'
                    ' calling genotype'),
    'author': 'Shujia Huang',
    'url': 'http://gitlab.biosoft.pub/huangshujia/BaseVar',
    'download_url': 'http://gitlab.biosoft.pub/huangshujia/BaseVar',
    'author_email': 'huangshujia9@gmail.com.',
    'version': '0.0.1.dev1',
    'install_requires': ['pysam'],
    'packages': ['basevar'],
    'scripts': [],
    'name': 'BaseVar'
}

setup(**config)
