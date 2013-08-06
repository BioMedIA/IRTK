#!/usr/bin/python

import irtk
import sys
import os
import numpy as np
import scipy.ndimage as nd
import cv2
import argparse

#from irtk.ext.patches import extract_oriented_patches2D as extract_patches2D
from irtk.ext.patches import extract_patches2D
from irtk.ext.crf import crf
from lib.BundledSIFT import get_OFD

from irtk.ext.opencv import sift_patches

from sklearn.externals import joblib
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(
    description='Slice-by-slice masking of fetal brain MRI (3D).' )
parser.add_argument( '--img', nargs='+', type=str, required=True )
parser.add_argument( '--mask', nargs='+', type=str, required=True )
parser.add_argument( '--ga', type=float, required=True )
parser.add_argument( '--output_dir', type=str, required=True )
parser.add_argument( '-r', '--radius', type=int, default=8,
                     help="patch size" )
parser.add_argument( '--cpu', type=int, default=-1,
                     help="number of CPUs used" )
parser.add_argument( '--ntrees', type=int, default=30,
                     help="number of trees" )
parser.add_argument( '--debug', action="store_true", default=False )
args = parser.parse_args()
print args

DEBUG = args.debug

if not os.path.exists( args.output_dir ):
    os.makedirs( args.output_dir )

def get_training_data( file_img, file_mask, r ):
    # create mask
    input_mask = irtk.imread( file_mask )
    x_min, y_min, z_min, x_max, y_max, z_max = (input_mask == 0).bbox()

    background = irtk.zeros( input_mask.get_header(), dtype='uint8' )
    background[z_min:z_max+1,
               y_min:y_max+1,
               x_min:x_max+1] = 1
    background = nd.morphological_gradient( background, size=7)
    n = background[z_min+1:z_max,
                   y_min+1:y_max,
                   x_min+1:x_max].sum()
    z = np.random.randint(low=0, high=input_mask.shape[0],size=1.25*n)
    y = np.random.randint(low=0, high=input_mask.shape[1],size=1.25*n)
    x = np.random.randint(low=0, high=input_mask.shape[2],size=1.25*n)
    background[z,y,x] = 1
    background[z_min+1:z_max,
               y_min+1:y_max,
               x_min+1:x_max] = 0
    
    foreground = (input_mask == 1).astype('uint8')

    new_mask = irtk.zeros( input_mask.get_header(), dtype='uint8' )
    new_mask[foreground == 1] = 1
    new_mask[background != 0] = 2

    img = irtk.imread( file_img, dtype='float32' )#.saturate().rescale().astype('uint8')

    X = []
    Y = []

    for z in xrange(img.shape[0]):
        YX = np.transpose( np.nonzero( foreground[z] ) )
        if DEBUG:
            YX = YX[::10]
        else:
            YX = YX[::2]
        if YX.shape[0] == 0:
            continue
        #patches = sift_patches( img[z], YX, 2*r+1 )
        patches = extract_patches2D( img[z], r, YX )
        patches = np.reshape( patches, (patches.shape[0],patches.shape[1]*patches.shape[2]) )
        print patches.shape, YX.shape
        #patches = np.hstack( (YX, [[z]]*len(YX), patches) )
        X.extend( patches )
        Y.extend( [1]*len(YX) )

    for z in xrange(img.shape[0]):
        YX = np.transpose( np.nonzero( background[z] ) )
        if DEBUG:
            YX = YX[::10]
        else:
            YX = YX[::2]
        if YX.shape[0] == 0:
            continue
        #patches = sift_patches( img[z], YX, 2*r+1 )
        patches = extract_patches2D( img[z], r, YX )
        patches = np.reshape( patches, (patches.shape[0],patches.shape[1]*patches.shape[2]) )
        print patches.shape, YX.shape
        #patches = np.hstack( (YX, [[z]]*len(YX), patches) )
        X.extend( patches )
        Y.extend( [0]*len(YX) )

    return X, Y

XY = Parallel(n_jobs=args.cpu)(delayed(get_training_data)(file_img, file_mask, args.radius)
                         for file_img, file_mask in zip(args.img,args.mask) )

print len(XY)

X = []
Y = []
for x,y in XY:
    X.extend(x)
    Y.extend(y)

X = np.array( X, dtype='float32').copy()
Y = np.array(Y, dtype='int32').copy()

n_positive = Y.sum()
n_points = Y.shape[0]
print "RATIO = ", n_positive, n_points, float(n_positive) / float(n_points) * 100
print "learning..."

from sklearn.ensemble import RandomForestClassifier
neigh = RandomForestClassifier( n_estimators=args.ntrees,
                                criterion='gini',
                                n_jobs=args.cpu )
neigh.fit(X, Y)
neigh.set_params(n_jobs=1)

def mask_image( file_img, file_mask, ga, r, neigh, output_dir ):
    img = irtk.imread( file_img, dtype='float32' )#.saturate().rescale().astype('uint8')
    input_mask = irtk.imread( file_mask )
    
    print "predicting..."
    res = irtk.zeros( img.get_header(), dtype='float32' )
    res2 = irtk.zeros( img.get_header(), dtype='float32' )
    res3 = irtk.zeros( img.get_header(), dtype='float32' )
    res4 = irtk.zeros( img.get_header(), dtype='uint8' )
    mask = irtk.ones( input_mask.get_header(), dtype='uint8' )
    mask[input_mask == 2] = 0
    for z in xrange(img.shape[0]):
        print z
        YX = np.transpose( np.nonzero( mask[z] ) )
        if YX.shape[0] == 0:
            continue # this slice does not intersect the box
        #YX = YX[::100]
        #patches = sift_patches( img[z], YX, 2*r+1 )
        patches = extract_patches2D( img[z], r, YX )
        patches = np.reshape( patches, (patches.shape[0],patches.shape[1]*patches.shape[2]) )
        #patches = np.hstack( (YX, [[z]]*len(YX), patches) ) 

        predictions = neigh.predict_proba(patches)[:,1]
        res[z,YX[:,0],YX[:,1]] = predictions

        x_min, y_min, x_max, y_max = mask[z].bbox()
        cropped_img = img[z,
                          y_min:y_max+1,
                          x_min:x_max+1].view(np.ndarray).astype("float32").copy()
        cropped_proba = res[z,
                            y_min:y_max+1,
                            x_min:x_max+1].view(np.ndarray).astype("float64").copy()
        cropped_labels = (cropped_proba > 0.5).astype("int32").copy()
    
        cropped_labels = crf( cropped_img,
                              cropped_labels,
                              cropped_proba,
                              l=1.0,
                          #degree=1
                              )
        res2[z,
             y_min:y_max+1,
             x_min:x_max+1] = cropped_labels

        # get connected component closest to center
        labeled_array, num_features = nd.label( res2[z] )
        
        if num_features == 0:
            print "no object found"
            continue        

        m = -1
        largest = -1
        for i in xrange( 1, num_features+1 ):
            s = np.sum( labeled_array == i )
            if s > m:
                m = s
                largest = i

        res3[z][labeled_array == largest] = 1

    # fill holes, close and dilate
    disk_close = irtk.disk( 5 )
    disk_dilate = irtk.disk( 2 )
    for z in xrange(mask.shape[0]):
        res3[z] = nd.binary_fill_holes( res3[z] )
        res3[z] = nd.gaussian_filter( res3[z], 4 ) > 0.5 
        #res3[z] = nd.binary_closing( res3[z], disk_close )
    #    res3[z] = nd.binary_dilation( res3[z], disk_dilate )
    
    # re-read the image in case we processed it
    img = irtk.imread( file_img, dtype='float32' )

    min_x, min_y, min_z, max_x, max_y, max_z = res3.bbox()
    x0 = min_x + (max_x - min_x) / 2
    y0 = min_y + (max_y - min_y) / 2
    ofd = get_OFD(ga)/img.header['pixelSize'][0]

    cropped_img = img[min_z:max_z+1,
                      max(0,int(round(y0-ofd/2))):min(img.shape[1],int(round(y0+ofd/2+1))),
                      max(0,int(round(x0-ofd/2))):min(img.shape[2],int(round(x0+ofd/2+1)))].copy()

    irtk.imwrite(output_dir + "/very_large_"+os.path.basename(file_img),
                 cropped_img )
    
    cropped_proba = res[min_z:max_z+1,
                        max(0,int(round(y0-ofd/2))):min(img.shape[1],int(round(y0+ofd/2+1))),
                        max(0,int(round(x0-ofd/2))):min(img.shape[2],int(round(x0+ofd/2+1)))].copy()

    irtk.imwrite(output_dir + "/proba_"+os.path.basename(file_img),
                 cropped_proba )

    img[res3 == 0] = -1
    masked = img.bbox(crop=True)
    irtk.imwrite(output_dir + "/masked_"+os.path.basename(file_img), masked )

XY = Parallel(n_jobs=args.cpu)(delayed(mask_image)(file_img, file_mask, args.ga,
                         args.radius, neigh, args.output_dir)
                         for file_img, file_mask in zip(args.img,args.mask) )
