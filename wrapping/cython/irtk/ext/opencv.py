import numpy as np
import _opencv

def sift_patches( img, YX, size ):
    img = img.astype('uint8').copy()
    YX = np.array( YX, dtype='int32' ).copy()
    return _opencv.sift_patches( img, YX, size )
