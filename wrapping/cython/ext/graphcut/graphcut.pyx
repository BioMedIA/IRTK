import cython
import numpy as np
cimport numpy as np

__all__ = ["graphcut"]

np.import_array()

ctypedef double voxel_t

ctypedef double captype
ctypedef double tcaptype
ctypedef double flowtype

cdef extern from "_graphcut.h":
    void _graphcut( voxel_t*,
                   int, int, int,
                   double,
                   unsigned char*,
                   double*,
                   unsigned char* )

@cython.boundscheck(False)
@cython.wraparound(False)     
def graphcut( np.ndarray[voxel_t, ndim=3, mode="c"] img,
              np.ndarray[unsigned char, ndim=3, mode="c"] mask,
              np.ndarray[double, ndim=1, mode="c"] raw_spacing ):

    cdef double img_std = np.std( img )

    cdef np.ndarray[unsigned char, ndim=3, mode="c"] seg = np.zeros( (img.shape[0],
                                                                      img.shape[1],
                                                                      img.shape[2]),
                                                                     dtype='uint8')

    print "starting graphcut..."                                                                 
    _graphcut( <voxel_t*> img.data,
               img.shape[0],
               img.shape[1],
               img.shape[2],
               img_std,
               <unsigned char*> mask.data,
               <double*> raw_spacing.data,
               <unsigned char*> seg.data )
    
    return seg
