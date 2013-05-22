from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy.distutils.misc_util import get_numpy_include_dirs

from glob import glob
import os
import sys

for i in xrange(1,len(sys.argv)):
    if sys.argv[i] == "--build-temp":
        build_tmp = sys.argv[i+1]
    if sys.argv[i] == "--build-lib":
        build_lib = sys.argv[i+1]

IRTK_src = "../../"
IRTK_build = build_lib + "/../../"

TMP_DIR = build_tmp

if os.environ['USER'] == "kevin":
    VTK_libdir = "/usr/lib/"
    VTK_include_dir = "/usr/include/vtk-5.8/"
    TBB_libdir = VTK_libdir
    TBB_include_dir = VTK_include_dir
    extra_link_args = ["-lpython2.7"]

else:
    VTK_libdir = "/vol/vipdata/users/kpk09/LOCAL/lib/vtk-5.8"
    VTK_include_dir = "/vol/vipdata/users/kpk09/LOCAL/include/vtk-5.8"
    TBB_libdir = "/vol/vipdata/users/kpk09/LOCAL/lib"
    TBB_include_dir = "/vol/vipdata/users/kpk09/LOCAL/include"
    extra_link_args = ["-lpython2.6"]

TBB_libs = ["tbb", "rt"]

# http://docs.python.org/2/distutils/apiref.html

# Concatenate automatically generated files
import fileinput
with open(TMP_DIR + "/tmp_irtk.pyx", 'w') as fout:
    for line in fileinput.input(["src/_irtk.pyx",TMP_DIR + "/templates.pyx"]):
        fout.write(line)
    fout.close()

def get_IRTK_include_dirs( folder=IRTK_src ):
    include_dirs = [ "recipes/include",
                     "common++/include",
                     "geometry++/include",
                     "image++/include",
                     "contrib++/include",
                     "recipes/include",
                     "packages/transformation/include",
                     "packages/registration/include",
                     "packages/registration2/include",
                     "packages/segmentation/include"]
    return map(
        lambda x: folder + "/" + x,
        include_dirs )

def get_IRTK_static_libraries( folder=IRTK_build ):
    static_libraries = [ "libznz.a",
                         "libniftiio.a",
                         "libcommon++.a",
                         "libsegmentation++.a",
                         "libimage++.a",
                         "librecipes.a",                         
                         "libtransformation++.a",                         
                         "libregistration++.a",
                         "libregistration2++.a",                         
                         "libcontrib++.a",
                         "libgeometry++.a",
                         "librview++.a" ]
    return map(
        lambda x: folder + "/lib/" + x,
        static_libraries )

def get_VTK_libs( folder=VTK_libdir ):
    libs = glob( folder + "/libvtk*.so" )
    return map( lambda x: os.path.basename( x )[3:-3], libs )
        
setup(
    cmdclass = { 'build_ext' : build_ext },
    ext_modules=[
        Extension( "_irtk",
                  [ 
                    "src/registration.cc",
                    "src/reconstruction.cc",
                    TMP_DIR + "/templates.cc",
                    "src/irtk2cython.cc",
                    "src/voxellise.cc",
                    "src/drawing.cc",
                    TMP_DIR + "/tmp_irtk.pyx"],
                   language="c++",
                   include_dirs = get_numpy_include_dirs()
                                + get_IRTK_include_dirs()
                                + [ "include", TMP_DIR, VTK_include_dir,
                                    TBB_include_dir ],
                   library_dirs = [ IRTK_build + "/lib", VTK_libdir, TBB_libdir ],
                   # the order of the libraries matters
                   libraries = [ "z","m","SM","dl","nsl","png","jpeg",
                                 "niftiio",
                                 "znz",                
                                 "common++",
                                 "segmentation++",
                                 "image++",
                                 "recipes",
                                 "transformation++",                                 
                                 "registration++",    
                                 "registration2++",
                                 "contrib++",
                                 "geometry++",
                                 "rview++" ]
                               + get_VTK_libs()
                               + TBB_libs, 
                   extra_objects = get_IRTK_static_libraries(),
                   extra_compile_args = ["-fPIC", "-O2", "-rdynamic",
                                         "-DHAS_TBB", "-DHAS_VTK", "-DHAS_VTK_HYBRID"],
                   runtime_library_dirs = [IRTK_build + "/lib"],
                   extra_link_args = [ "-Wl,-no-undefined", "-rdynamic" ]
                                     + extra_link_args
                   ),
        
        # ext directory (non IRTK stuff)
        Extension( "_template",
                   [ "ext/_template.pyx" ],
                   #language="c++",
                   include_dirs = get_numpy_include_dirs()
                   ),

        Extension( "_slic",
                   [ "ext/SLIC/_slic.pyx", "ext/SLIC/SLIC.cc" ],
                   language = "c++",
                   include_dirs = get_numpy_include_dirs()
                                + [ "ext/SLIC" ]
                   ),

        Extension("graphcut", [ "ext/graphcut/graphcut.pyx",
                                "ext/graphcut/_graphcut.cc",
                                "ext/graphcut/maxflow-v3.02.src/graph.cpp",
                                "ext/graphcut/maxflow-v3.02.src/maxflow.cpp" ],
                  language="c++",
                  include_dirs = get_numpy_include_dirs()
                               + [ "ext/graphcut/",
                                   "ext/graphcut/maxflow-v3.02.src" ]
                  ),
        ]
    )


