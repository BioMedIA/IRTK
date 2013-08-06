__all__ = [ "voxellise",
            "shrinkDisk" ]

import image
import _irtk
import numpy as np

def voxellise( points, triangles, header=None, pixelSize=[1,1,1,1] ):
    points = np.array( points, dtype='float64' )
    triangles = np.array( triangles, dtype='int32' )
    
    if header is None:
        pixelSize = np.array( pixelSize, dtype='float64')
        x_min, y_min, z_min = points.min(axis=0) - pixelSize[:3]
        x_max, y_max, z_max = points.max(axis=0) + pixelSize[:3]
        origin = [x_min + (x_max-x_min)/2,
                  y_min + (y_max-y_min)/2,
                  z_min + (z_max-z_min)/2,
                  0]
        dim = [ (x_max - x_min) / pixelSize[0] + 1,
                (y_max - y_min) / pixelSize[1] + 1,
                (z_max - z_min) / pixelSize[2] + 1,
                1]
        header = image.new_header( origin=origin,
                                   dim=dim,
                                   pixelSize=pixelSize )

    img = _irtk.voxellise( points, triangles, header )

    return image.Image( img, header )


def shrinkDisk( img,
                center=None,
                radius=None,
                steps=50 ):
    if center is None:
        center = np.array([float(img.shape[0])/2,
                           float(img.shape[1])/2],
                          dtype='float64')
    if radius is None:
        radius = float(img.shape[0])/2
    img = img.astype('uint8').copy()
    center = np.array( center, dtype='float64' ).copy()
    return _irtk.shrinkDisk(img, center, radius, steps )
    
