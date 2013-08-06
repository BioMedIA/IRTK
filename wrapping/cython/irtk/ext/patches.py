import numpy as np
import _patches
import math
import scipy.ndimage as nd
from numbers import Number
import cv2

def pad2D(x,padding=(1,1,0,0)):
    if len(x.shape) == 2:
        padded = np.zeros( (x.shape[0]+padding[0]+padding[2],
                            x.shape[1]+padding[1]+padding[3]),
                           dtype='float64' )
    elif len(x.shape) == 3:
        padded = np.zeros( (x.shape[0]+padding[0]+padding[2],
                            x.shape[1]+padding[1]+padding[3],
                            x.shape[2]),
                           dtype='float64' )
    else:
        raise ValueError("wrong number of dimensions")
    padded[padding[0]:padding[0]+x.shape[0],
           padding[1]:padding[1]+x.shape[1]] = x
    return padded

def integral_image2D(x):
    return pad2D( x.cumsum(1).cumsum(0) )

def extract_patches2D( img, size, coordinates ):
    if isinstance(size, Number):
        size = (size,size)
    (rx, ry) = size
    X = coordinates[:,1].astype('int32').copy()
    Y = coordinates[:,0].astype('int32').copy()
    return _patches.extract_patches2D( img.astype('float32').copy(),
                                       X, Y,
                                       rx, ry )

def reconstruct_from_patches2D( shape, patches, yx, average=True ):
    img = np.zeros( shape, dtype="float32")
    rx = (patches.shape[2]-1) / 2
    ry = (patches.shape[1]-1) / 2
    if average:
        for i in xrange(patches.shape[0]):
            img[yx[i,0]-ry:yx[i,0]+ry+1,
                yx[i,1]-rx:yx[i,1]+rx+1] += patches[i]
    else:
        for i in xrange(patches.shape[0]):
            img[yx[i,0]-ry:yx[i,0]+ry+1,
                yx[i,1]-rx:yx[i,1]+rx+1] = patches[i]
    return img

def extract_oriented_patches2D( img, r, coordinates, nb_proj=100 ):

    img = img.astype('float64')

    projection_matrix = np.zeros( (nb_proj,2),
                                  dtype='float64' )
    for i in xrange(nb_proj):
        theta = float(i) * 2.0 * math.pi / nb_proj
        projection_matrix[i,0] = math.sin(theta)
        projection_matrix[i,1] = math.cos(theta)

    print "computing gradients..."
    # grad = np.dstack( ( nd.sobel( img, mode='constant', axis=0 ),
    #                     nd.sobel( img, mode='constant', axis=1 ) ) )

    grad = np.dstack( ( nd.convolve1d( img, [-1,1], mode='constant', axis=0 ),
                        nd.convolve1d( img, [-1,1], mode='constant', axis=1 ) ) )
    
    print "projecting gradients..."
    hist = _patches.project_gradients2D( grad, projection_matrix )
    hist = integral_image2D( hist )

    print hist

    print "extracting patches..."
    Y = coordinates[:,0].copy().astype('int32')
    X = coordinates[:,1].copy().astype('int32')
    return _patches.extract_oriented_patches2D( img,
                                                hist,
                                                projection_matrix,
                                                X,
                                                Y,
                                                r )

def sift_patches( img, size, YX ):
    keypoints = []
    for y,x in YX:
        keypoints.append( cv2.KeyPoint( x, y, size ) )
    siftExtractor = cv2.DescriptorExtractor_create("SIFT")
    (keypoints, descriptors) = siftExtractor.compute(img,keypoints)
    return descriptors
