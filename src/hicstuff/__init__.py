"""Stuff to do with Hi-C data

A simple library/pipeline to generate and handle simple Hi-C data.

"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("hicstuff")
except PackageNotFoundError:
    __version__ = "unknown"

from .hicstuff import *

__author__ = "Lyam Baudry, \
              Axel Cournac, \
              Cyril Matthey-Doret, \
              Rémi Montagne, \
              Vittore Scolari"
__copyright__ = "Copyright © 2017-2019, Institut Pasteur, Paris, France"
__credits__ = [
    "Lyam Baudry",
    "Axel Cournac",
    "Cyril Matthey-Doret",
    "Rémi Montagne",
    "Vittore Scolari",
]
__license__ = "GPLv3"
__maintainer__ = "Lyam Baudry, \
                  Cyril Matthey-Doret"
__email__ = "lyam.baudry@pasteur.fr, cmdoret@pasteur.fr"
__status__ = "Beta"
