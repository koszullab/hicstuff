#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""General purpose stuff to generate and handle Hi-C data in its simplest form.
"""

from setuptools import setup, find_packages
import codecs

CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

name = "hicstuff"

MAJOR = 3
MINOR = 1
MAINTENANCE = 1
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MAINTENANCE)

LICENSE = "GPLv3"
URL = "https://github.com/koszullab/hicstuff"

DESCRIPTION = __doc__.strip("\n")

with codecs.open("README.md", encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

with open("requirements.txt", "r") as f:
    REQUIREMENTS = f.read().splitlines()

with open("hicstuff/version.py", "w") as f:
    f.write("__version__ = '{}'\n".format(VERSION))


setup(
    name=name,
    author="cyril.matthey-doret@pasteur.fr",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    version=VERSION,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    url=URL,
    packages=find_packages(),
    # package_data={"hicstuff": ("kernels/*")},
    python_requires=">=3.6",
    include_package_data=True,
    install_requires=REQUIREMENTS,
    long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["hicstuff=hicstuff.main:main"]},
    extras_require={"mappy": ["mappy"], "cooler": ["cooler"]},
)
