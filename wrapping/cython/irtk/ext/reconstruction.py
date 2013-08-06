__all__ = [ "reconstruct" ]

import numpy as np
import _irtk

import image
from registration import RigidTransformation

def reconstruct( slice_list,
                 transformation_list,
                 stack_ids,
                 mask,
                 thickness,
                 iterations=1 ):
    """
    No registration is performed if there is only one iteration.
    
    """

    # if mask is None:
    #     mask = image.ones( template_header, dtype='float32' )
        
    stack_ids = np.array( stack_ids, dtype='int32' )
    thickness = np.array( thickness, dtype='float64' )
    
    n = len(slice_list)
    nb_pixels = 0
    for img in slice_list:
        nb_pixels += img.nb_pixels()

    pixelData = np.zeros( nb_pixels, dtype='float32' )
    pixelSize = np.zeros( 4*n, dtype='float64' )
    xAxis = np.zeros( 3*n, dtype='float64' )
    yAxis = np.zeros( 3*n, dtype='float64' )
    zAxis = np.zeros( 3*n, dtype='float64' )
    origin = np.zeros( 4*n, dtype='float64' )
    dim = np.zeros( 4*n, dtype='int32' )

    tx = np.zeros( n, dtype='float64' )
    ty = np.zeros( n, dtype='float64' )
    tz = np.zeros( n, dtype='float64' )
    rx = np.zeros( n, dtype='float64' )
    ry = np.zeros( n, dtype='float64' )
    rz = np.zeros( n, dtype='float64' )
    
    offset = 0
    for i, img in enumerate(slice_list):
        pixelData[offset:offset+img.nb_pixels()] = img.flatten()
        offset += img.nb_pixels()
        pixelSize[4*i:4*(i+1)] = img.header['pixelSize']
        xAxis[3*i:3*(i+1)] = img.header['orientation'][0]
        yAxis[3*i:3*(i+1)] = img.header['orientation'][1]
        zAxis[3*i:3*(i+1)] = img.header['orientation'][2]
        origin[4*i:4*(i+1)] = img.header['origin']
        dim[4*i:4*(i+1)] = img.header['dim']

        tx[i] = transformation_list[i].tx
        ty[i] = transformation_list[i].ty
        tz[i] = transformation_list[i].tz
        rx[i] = transformation_list[i].rx
        ry[i] = transformation_list[i].ry
        rz[i] = transformation_list[i].rz

    #_irtk.write_list( pixelData, pixelSize, xAxis, yAxis, zAxis, origin, dim, n )        

    reconstructed = _irtk.reconstruct( pixelData,
                                       pixelSize,
                                       xAxis,
                                       yAxis,
                                       zAxis,
                                       origin,
                                       dim,
                                       n,
                                       stack_ids,
                                       iterations,
                                       tx, ty, tz,
                                       rx, ry, rz,
                                       thickness,
                                       mask.get_data('float32','cython'),
                                       mask.get_header() )

    new_transformations = []
    for i in xrange(n):
        new_transformations.append( RigidTransformation( tx=tx[i], ty=ty[i], tz=tz[i],
                                                         rx=rx[i], ry=ry[i], rz=rz[i] ) )

    return ( image.Image(reconstructed, mask.get_header()), new_transformations )
