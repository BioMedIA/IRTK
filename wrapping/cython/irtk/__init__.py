"""
This module is a wrapper around IRTK (http://www.doc.ic.ac.uk/~dr/software/)
written using Cython and a script to simulate templated code.
"""

import sys

__all__ = []

import irtk.image
from irtk.image import *
import irtk.registration
from irtk.registration import *
import irtk.cmd
from irtk.cmd import *
import irtk.morphology
from irtk.morphology import *
__all__.extend(irtk.image.__all__)
__all__.extend(irtk.registration.__all__)
__all__.extend(irtk.cmd.__all__)
__all__.extend(irtk.morphology.__all__)
