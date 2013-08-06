from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy.distutils.misc_util import get_numpy_include_dirs
        
setup(
    cmdclass = { 'build_ext' : build_ext },
    ext_modules=[
       
        # ext directory (non IRTK stuff)
        Extension( "_template",
                   [ "_template.pyx" ],
                   #language="c++",
                   include_dirs = get_numpy_include_dirs()
                   ),

        Extension( "_slic",
                   [ "SLIC/_slic.pyx", "SLIC/SLIC.cc" ],
                   language = "c++",
                   include_dirs = get_numpy_include_dirs()
                                + [ "SLIC" ]
                   ),

        Extension("graphcut", [ "graphcut/graphcut.pyx",
                                "graphcut/_graphcut.cc",
                                "graphcut/maxflow-v3.02.src/graph.cpp",
                                "graphcut/maxflow-v3.02.src/maxflow.cpp" ],
                  language="c++",
                  include_dirs = get_numpy_include_dirs()
                               + [ "graphcut/",
                                   "graphcut/maxflow-v3.02.src" ]
                  ),

        Extension("_clahe", [ "CLAHE/_clahe.pyx",
                             "CLAHE/clahe.c" ],
                  language="c",
                  include_dirs = get_numpy_include_dirs()
                               + [ "CLAHE/" ]
                  ),

        Extension("crf", [ "crf/crf.pyx",
                           "crf/_crf.cc",
                           "crf/src/graph.cpp",
                           "crf/src/GCoptimization.cpp",
                           "crf/src/LinkedBlockList.cpp",
                           "crf/src/maxflow.cpp" ],
                  language="c++",
                  include_dirs = get_numpy_include_dirs()
                               + [ "crf/",
                                   "crf/src" ]
                  ),

        Extension("_patches", [ "patches/_patches.pyx",
                               "patches/patches.cc" ],
                  language = "c++",
                  include_dirs = get_numpy_include_dirs()
                                 + [ "patches/" ]
                  ),
        ]
    )


