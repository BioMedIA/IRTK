import sys
if "matplotlib" not in sys.modules.keys():
    # then we do not have an interactive pylab environment
    import matplotlib
    matplotlib.use('PDF') # pygtk not installed on biomedia machines

__all__ = []

import _irtk
_irtk.initialise_library()

# core
import irtk.image
from irtk.image import *
import irtk.registration
from irtk.registration import *
import irtk.vtk2irtk
from irtk.vtk2irtk import *
import irtk.reconstruction
from irtk.reconstruction import *
import irtk.cmd
from irtk.cmd import *
__all__.extend(irtk.image.__all__)
__all__.extend(irtk.registration.__all__)
__all__.extend(irtk.vtk2irtk.__all__)
__all__.extend(irtk.reconstruction.__all__)
__all__.extend(irtk.cmd.__all__)

# ext
import irtk.ext.template
from irtk.ext.template import *
import irtk.ext.Show3D
from irtk.ext.Show3D import *
import irtk.ext.slic
from irtk.ext.slic import *
import irtk.ext.graphcut
from irtk.ext.graphcut import *
__all__.extend(irtk.ext.template.__all__)
__all__.extend(irtk.ext.Show3D.__all__)
__all__.extend(irtk.ext.slic.__all__)
__all__.extend(irtk.ext.graphcut.__all__)
