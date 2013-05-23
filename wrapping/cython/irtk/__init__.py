import sys
if "matplotlib" not in sys.modules.keys():
    # then we do not have an interactive pylab environment
    import matplotlib
    matplotlib.use('PDF') # pygtk not installed on biomedia machines

try:
    import cv2 # OpenCV
    has_opencv = True
except ImportError:
    has_opencv = False

try:
    import vtk # VTK
    has_vtk = True
except ImportError:
    has_vtk = False    

# http://stackoverflow.com/questions/7749914/add-method-to-loaded-class-module-in-python

__all__ = []

import _irtk
_irtk.initialise_library()

# core
import irtk.image
from irtk.image import *
import irtk.registration
from irtk.registration import *
import irtk.reconstruction
from irtk.reconstruction import *
import irtk.cmd
from irtk.cmd import *
__all__.extend(irtk.image.__all__)
__all__.extend(irtk.registration.__all__)
__all__.extend(irtk.reconstruction.__all__)
__all__.extend(irtk.cmd.__all__)

# ext
import irtk.ext.template
from irtk.ext.template import *
import irtk.ext.slic
from irtk.ext.slic import *
import irtk.ext.graphcut
from irtk.ext.graphcut import *
__all__.extend(irtk.ext.template.__all__)
__all__.extend(irtk.ext.slic.__all__)
__all__.extend(irtk.ext.graphcut.__all__)

if has_vtk:
    import irtk.vtk2irtk
    from irtk.vtk2irtk import *
    __all__.extend(irtk.vtk2irtk.__all__)
    import irtk.ext.Show3D
    from irtk.ext.Show3D import *
    __all__.extend(irtk.ext.Show3D.__all__)

if has_opencv:
    import irtk.opencv_functions
    from irtk.opencv_functions import *
    __all__.extend(irtk.opencv_functions.__all__)
    irtk.Image.resample2D = irtk.opencv_functions.resample2D
