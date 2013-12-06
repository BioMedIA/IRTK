#!/usr/bin/python

import irtk
import sys
import os
import numpy as np
import scipy.ndimage as nd
import cv2
import argparse

#from skimage.filter import denoise_tv_chambolle 

#from irtk.ext.patches import extract_oriented_patches2D as extract_patches2D
from irtk.ext.patches import extract_patches2D
from lib.BundledSIFT import get_OFD

from sklearn.externals import joblib
from joblib import Parallel, delayed


from scipy.stats.mstats import mquantiles
    
parser = argparse.ArgumentParser(
    description='Slice-by-slice masking of fetal brain MRI (3D).' )
parser.add_argument( '--img', nargs='+', type=str, required=True )
parser.add_argument( '--mask', nargs='+', type=str, required=True )
parser.add_argument( '--ga', type=float, required=True )
parser.add_argument( '--output_dir', type=str, required=True )
parser.add_argument( '-r', '--radius', type=int, default=8,
                     help="patch size" )
parser.add_argument( '-l', '--l', type=float, default=30.0,
                     help="lambda" )
parser.add_argument( '--cpu', type=int, default=-1,
                     help="number of CPUs used" )
parser.add_argument( '--ntrees', type=int, default=30,
                     help="number of trees" )
parser.add_argument( '--do_3D', action="store_true", default=False )
parser.add_argument( '--do_patchZ', action="store_true", default=False )
parser.add_argument( '--no_cleaning', action="store_true", default=False )
parser.add_argument( '--debug', action="store_true", default=False )
parser.add_argument( '--mass', action="store_true", default=False )
args = parser.parse_args()
print args

DEBUG = args.debug

if not os.path.exists( args.output_dir ):
    os.makedirs( args.output_dir )

def get_BV( GA ):
    """
    Return expected brain volume according to gestational age.

    Reference:
    "The assessment of normal fetal brain volume by 3-D ultrasound"
    Chiung-Hsin Chang, Chen-Hsiang Yu, Fong-Ming Chang, Huei-Chen Ko, Hsi-Yao Chen
    """
    # mL to mm3 , 1 ml is 1 cm^3
    return (-171.48036 + 4.8079*GA + 0.29521*GA**2)*1000

def get_noiseXY(img):
    img = img.astype('float32')
    new_img = np.zeros(img.shape,dtype='float32')
    for z in xrange(img.shape[0]):
        new_img[z] = nd.gaussian_filter( img[z], 2, mode='reflect' )
    noise = img - new_img
    #print "Noise XY:", noise.std(), img.std()
    return noise.std()

def get_noiseZ(img):
    img = img.astype('float32')
    new_img = np.zeros(img.shape,dtype='float32')
    for x in xrange(img.shape[2]):
        new_img[:,:,x] = nd.gaussian_filter( img[:,:,x], 2, mode='reflect' )
    noise = img - new_img
    #print "Noise Z:", noise.std(), img.std()
    return noise.std()

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

    img = irtk.imread( file_img, dtype='float32' )
    
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
        patches = extract_patches2D( img[z], r, YX )
        patches = np.reshape( patches, (patches.shape[0],patches.shape[1]*patches.shape[2]) )
        print patches.shape, YX.shape
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
        patches = extract_patches2D( img[z], r, YX )
        patches = np.reshape( patches, (patches.shape[0],patches.shape[1]*patches.shape[2]) )
        print patches.shape, YX.shape
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
    img = irtk.imread( file_img, dtype='float32' )

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
        patches = extract_patches2D( img[z], r, YX )
        patches = np.reshape( patches, (patches.shape[0],patches.shape[1]*patches.shape[2]) )

        predictions = neigh.predict_proba(patches)[:,1]
        res[z,YX[:,0],YX[:,1]] = predictions

    x_min, y_min, z_min, x_max, y_max, z_max = mask.bbox()

    proba = res[z_min:z_max+1,
                y_min:y_max+1,
                x_min:x_max+1]

    if args.mass:
        BV = get_BV( args.ga )
        box_volume = (z_max-z_min)*img.header['pixelSize'][2]*(y_max-y_min)*img.header['pixelSize'][1]*(x_max-x_min)*img.header['pixelSize'][0]
        ratio = float(BV) / float(box_volume)
        print "ratio", ratio
        q0,q1 = mquantiles( proba.flatten(), prob=[0.5*(1.0-ratio),
                                                   1.0-0.5*ratio] )
        print "threshold", q0,q1
        #threshold = max(0.5,threshold)
    
        # labels = res[z_min:z_max+1,
        #              y_min:y_max+1,
        #              x_min:x_max+1] > threshold
        
    #res = 1 / (np.exp(-(res-threshold)/(res.max()-res.min())))

        res[res<q0] = q0
        res[res>q1] = q1
        res -= res.min()
        res /= res.max()

    labels = res[z_min:z_max+1,
                 y_min:y_max+1,
                 x_min:x_max+1] > 0.5
   
    proba = res[z_min:z_max+1,
                y_min:y_max+1,
                x_min:x_max+1]
    
    cropped_img = img[z_min:z_max+1,
                      y_min:y_max+1,
                      x_min:x_max+1]

    if args.do_3D:
        labels = irtk.crf( cropped_img,
                           labels,
                           proba,
                           l=args.l,
                           sigma=get_noiseXY(cropped_img),
                           sigmaZ=get_noiseZ(cropped_img) )
    # elif args.do_patchZ:
    #     labels = irtk.crf_patchZ( cropped_img,
    #                               labels,
    #                               proba,
    #                               l=10.0 )   
    # else:
    #     for z in xrange(z_min,z_max+1):
    #         labels[z] = irtk.crf( cropped_img[z],
    #                               labels[z],
    #                               proba[z],
    #                               l=1.0 )

    print "MAX LABEL:", labels.max()
    irtk.imwrite(output_dir + "/bare_"+os.path.basename(file_img), labels )
    tmp = irtk.zeros( img.get_header(), dtype='uint8' )
    tmp[z_min:z_max+1,
        y_min:y_max+1,
        x_min:x_max+1] = labels
    ( min_x_bare, min_y_bare, min_z_bare,
      max_x_bare, max_y_bare, max_z_bare ) = tmp.bbox()
    
    if not args.no_cleaning:
        # clean by fitting ellipses enlarged of 10%
        for z in xrange(labels.shape[0]):
            edges = nd.morphological_gradient( labels[z] > 0,size=5 )
            points = np.transpose(edges.nonzero())[:,::-1]
            if len(points) == 0:
                continue
            points = np.array(map(lambda x:[x],points),dtype='int32')
            ellipse = cv2.fitEllipse(points)
            cv2.ellipse( labels[z], (ellipse[0],
                                     (1.1*ellipse[1][0],1.1*ellipse[1][1]),
                                     ellipse[2]) , 1, -1 )

    irtk.imwrite(output_dir + "/seg_"+os.path.basename(file_img), labels )
    irtk.imwrite(output_dir + "/res_"+os.path.basename(file_img), res )

    # re-read the image in case we processed it
    img = irtk.imread( file_img, dtype='float32' )
    cropped_img = img[z_min:z_max+1,
                      y_min:y_max+1,
                      x_min:x_max+1]
    cropped_img[labels==0] = -1
    masked = cropped_img.bbox(crop=True)
    irtk.imwrite(output_dir + "/masked_"+os.path.basename(file_img), masked )

    # re-read the image in case we processed it
    img = irtk.imread( file_img, dtype='float32' )    
    x0 = min_x_bare + (max_x_bare - min_x_bare) / 2
    y0 = min_y_bare + (max_y_bare - min_y_bare) / 2
    ofd = get_OFD(ga)/img.header['pixelSize'][0]

    cropped_img = img[min_z_bare:max_z_bare+1,
                      max(0,int(round(y0-ofd/2))):min(img.shape[1],int(round(y0+ofd/2+1))),
                      max(0,int(round(x0-ofd/2))):min(img.shape[2],int(round(x0+ofd/2+1)))].copy()

    irtk.imwrite(output_dir + "/very_large_"+os.path.basename(file_img),
                 cropped_img )
    
    cropped_proba = res[min_z_bare:max_z_bare+1,
                        max(0,int(round(y0-ofd/2))):min(img.shape[1],int(round(y0+ofd/2+1))),
                        max(0,int(round(x0-ofd/2))):min(img.shape[2],int(round(x0+ofd/2+1)))].copy()

    irtk.imwrite(output_dir + "/proba_"+os.path.basename(file_img),
                 cropped_proba )    
    

XY = Parallel(n_jobs=args.cpu)(delayed(mask_image)(file_img, file_mask, args.ga,
                         args.radius, neigh, args.output_dir)
                         for file_img, file_mask in zip(args.img,args.mask) )
