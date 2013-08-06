import cython
import numpy as np
cimport numpy as np

__all__ = ["crf"]

np.import_array()

ctypedef float pixel_t

ctypedef int LabelID
ctypedef float EnergyType

cdef extern from "_crf.h":
    void _crf( pixel_t* img,
               int shape0, int shape1,
               double std,
               LabelID* labels,
               double* proba,
               double l,
               int degree )
    void _crf3D( pixel_t* img,
          int shape0, int shape1, int shape2,
          double std,
          LabelID* labels,
          double* proba,
          double l )

@cython.boundscheck(False)
@cython.wraparound(False)     
def crf( np.ndarray[pixel_t, ndim=2, mode="c"] img,
         np.ndarray[LabelID, ndim=2, mode="c"] labels,
         np.ndarray[double, ndim=2, mode="c"] proba,
         double l=1.0,
         int degree=1 ):
    cdef double std = np.std( img )

    print "starting crf..."                                                                 
    _crf( <pixel_t*> img.data,
           img.shape[0],
           img.shape[1],
           std,
           <LabelID*> labels.data,
           <double*> proba.data,
           l,
           degree )
    
    return labels

@cython.boundscheck(False)
@cython.wraparound(False)     
def crf3D( np.ndarray[pixel_t, ndim=3, mode="c"] img,
         np.ndarray[LabelID, ndim=3, mode="c"] labels,
         np.ndarray[double, ndim=3, mode="c"] proba,
         double l=1.0 ):
    cdef double std = np.std( img )

    print "starting crf..."                                                                 
    _crf3D( <pixel_t*> img.data,
           img.shape[0],
           img.shape[1],
           img.shape[2],
           std,
           <LabelID*> labels.data,
           <double*> proba.data,
           l )
    
    return labels
