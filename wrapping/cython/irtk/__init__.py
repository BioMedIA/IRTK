import sys

if "matplotlib" not in sys.modules.keys():
    # then we do not have an interactive pylab environment
    import matplotlib
    matplotlib.use('PDF') # pygtk not installed on biomedia machines 

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
