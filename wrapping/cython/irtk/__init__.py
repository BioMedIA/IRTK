import sys
if "matplotlib" not in sys.modules.keys():
    # then we do not have an interactive pylab environment
    import matplotlib
    matplotlib.use('PDF') # pygtk not installed on biomedia machines

from irtk.Image import Image, imread, imwrite, zeros, ones, rview, display, imshow
from irtk.Transformation import RigidTransformation
from irtk.ext.template import match_template
#from irtk.ext.Show import imshow
from irtk.ext.Show3D import render
from irtk.ext.slic import slic

__all__ = [ "Image", "imread", "imwrite", "zeros", "ones", "rview", "display",
            "imshow",

            "RigidTransformation",
            
            "match_template",
            
            "render",

            "slic" ]
