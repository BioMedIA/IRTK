import cython
import numpy as np
cimport numpy as np

np.import_array()

ctypedef unsigned char kz_pixel_t

cdef extern from "clahe.h":
    int _CLAHE( kz_pixel_t* pImage,
                unsigned int uiXRes,
                unsigned int uiYRes,
                kz_pixel_t Min,
                kz_pixel_t Max,
                unsigned int uiNrX,
                unsigned int uiNrY,
                unsigned int uiNrBins,
                float fCliplimit)

def CLAHE(  np.ndarray[kz_pixel_t, ndim=2, mode="c"] img,
            unsigned int nrx,
            unsigned int nry,
            unsigned int nbins,
            float clip ):
    """
    pImage - Pointer to the input/output image
    uiXRes - Image resolution in the X direction
    uiYRes - Image resolution in the Y direction
    Min - Minimum greyvalue of input image (also becomes minimum of output image)
    Max - Maximum greyvalue of input image (also becomes maximum of output image)
    uiNrX - Number of contextial regions in the X direction (min 2, max uiMAX_REG_X)
    uiNrY - Number of contextial regions in the Y direction (min 2, max uiMAX_REG_Y)
    uiNrBins - Number of greybins for histogram ("dynamic range")
    float fCliplimit - Normalized cliplimit (higher values give more contrast)
  The number of "effective" greylevels in the output image is set by uiNrBins; selecting
  a small value (eg. 128) speeds up processing and still produce an output image of
  good quality. The output image will have the same minimum and maximum value as the input
  image. A clip limit smaller than 1 results in standard (non-contrast limited) AHE.

  Sources:
    http://tog.acm.org/resources/GraphicsGems/gemsiv/clahe.c
    http://public.cranfield.ac.uk/c5354/teaching/dip/opencv/lecture_demos/c/
    https://github.com/joshdoe/opencv-clahe/
    
    """
    cdef unsigned int uiXRes = img.shape[1]
    cdef unsigned int uiYRes= img.shape[0]
    cdef kz_pixel_t Min = np.min( img )
    cdef kz_pixel_t Max = np.max( img )

    error_messages = {
        -1 : "# of regions x-direction too large",
         -2 : "# of regions y-direction too large",
         -3 : "x-resolution no multiple of uiNrX",
         -4 : "y-resolution no multiple of uiNrY",
         -5 : "maximum too large",
         -6 : "minimum equal or larger than maximum",
         -7 : "at least 4 contextual regions required",
         -8 : "Not enough memory! (try reducing uiNrBins)"
         }

    cdef int error = _CLAHE( <kz_pixel_t*> img.data,
                              uiXRes,
                              uiYRes,
                              Min,
                              Max,
                              nrx,
                              nry,
                              nbins,
                              clip)

    if error is not 0:
        print error_messages[error]
        return None
    else:
        return img
