#!/usr/bin/python

import cv2
import sys,os
import csv
import numpy as np
from math import cos,sin
import math
import irtk
import argparse
import ransac

from glob import glob
from sklearn import neighbors
from sklearn import svm
from sklearn.externals import joblib
from joblib import Parallel, delayed

######################################################################

def is_in_ellipse( (x,y), ((xe,ye),(we,he),theta)):
    theta = theta / 180 * np.pi
    u = cos(theta)*(x-xe)+ sin(theta)*(y-ye)
    v = -sin(theta)*(x-xe)+cos(theta)*(y-ye)

    # http://answers.opencv.org/question/497/extract-a-rotatedrect-area/
    # http://felix.abecassis.me/2011/10/opencv-rotation-deskewing/
    # if theta < -45:
    #     tmp = we
    #     we = he
    #     he = we

    a = we/2
    b = he/2

    return (u/a)**2 + (v/b)**2 <= 1

def get_OFD( ga ):
    ofd_model = np.array([ -4.97315445e-03,
                            3.19846853e-01,
                            -2.60839214e+00,
                            2.62679565e+01])
    OFD = 0
    for k in range(4):
        OFD += ofd_model[3-k]*ga**k
    return OFD
    
def get_BPD( ga ):
    bpd_model = np.array([ -3.21975764e-03,
                            2.17264994e-01,
                            -1.81364324e+00,
                            2.00444225e+01])
    BPD = 0
    for k in range(4):
        BPD += bpd_model[3-k]*ga**k
    return BPD

def detect_mser( raw_file,
                 ga,
                 vocabulary,
                 mser_detector,
                 NEW_SAMPLING,
                 output_folder="debug",
                 DEBUG=False,
                 return_image_regions=False ):

    OFD = get_OFD(ga) / NEW_SAMPLING
    BPD = get_BPD(ga) / NEW_SAMPLING

    max_e = 0.64
    
    mser = cv2.MSER( _delta=5,
                     _min_area=60,
                     _max_area=14400,
                     _max_variation=0.15,
                     _min_diversity=.1,
                     _max_evolution=200,
                     _area_threshold=1.01,
                     _min_margin=0.003,
                     _edge_blur_size=5)

    sift = cv2.SIFT( nfeatures=0,
                     nOctaveLayers=3,
                     contrastThreshold=0.04,
                     edgeThreshold=10,
                     sigma=0.8)
    siftExtractor = cv2.DescriptorExtractor_create("SIFT")

    voca = np.load( open(vocabulary, 'rb') )
    classifier = neighbors.NearestNeighbors(1,algorithm='kd_tree')
    N = voca.shape[0] 
    classifier.fit(voca)
    #flann = pyflann.FLANN()
    #flann.build_index( voca.astype('float32') )


    svc = joblib.load(mser_detector)
    
    img = irtk.imread( raw_file, dtype='float32' )
    img = img.resample2D( NEW_SAMPLING ).saturate().rescale().astype('uint8')

    detections = []
    image_regions = []
    for z in range(img.shape[0]):
        detections.append([])
        image_regions.append([])
        
        # Extract MSER
        #print "extracting mser"
        contours = mser.detect(img[z,:,:])
        #print "mser done"
        
        if DEBUG:
            img_color = cv2.cvtColor( img[z,:,:], cv2.cv.CV_GRAY2RGB )
            for c in contours:
                ellipse = cv2.fitEllipse(np.array(map(lambda x:[x],
                                                  c),dtype='int32'))
                cv2.ellipse( img_color, (ellipse[0],
                                         (ellipse[1][0],ellipse[1][1]),
                                         ellipse[2]) , (0,0,255))

            cv2.imwrite(output_folder + "/" +str(z) + "_all_mser_.png",img_color )

            img_color_mser = cv2.cvtColor( img[z,:,:], cv2.cv.CV_GRAY2RGB )

        # Filter MSER
        selected_mser = []
        mask = np.zeros( (img.shape[1],img.shape[2]), dtype='uint8' )
        #print "fitting ellipses"
        for c in contours:
            ellipse = cv2.fitEllipse(np.reshape(c, (c.shape[0],1,2) ).astype('int32'))

            # filter by size
            if ( ellipse[1][0] > OFD
                 or ellipse[1][1] > OFD
                 or ellipse[1][0] < 0.5*OFD
                 or ellipse[1][1] < 0.5*OFD ) :
                continue
            
            # filter by eccentricity
            if math.sqrt(1-(np.min(ellipse[1])/np.max(ellipse[1]))**2) > max_e:
                continue

            cv2.ellipse( mask, ellipse, 255, -1 )
            selected_mser.append((c,ellipse))

        #print "ellipses done"
        if len(selected_mser) == 0:
            continue

        # Extract SIFT
        #print "extracting SIFT"
        keypoints = sift.detect(img[z,:,:],mask=mask)
        #print "SIFT done"
        if keypoints is None or len(keypoints) == 0:
            continue
        (keypoints, descriptors) = siftExtractor.compute(img[z,:,:],keypoints)

        # words = np.zeros(len(keypoints),dtype="int")
        # for i,d in enumerate(descriptors):
        #     words[i] = classifier.kneighbors(d, return_distance=False)
        words = classifier.kneighbors(descriptors, return_distance=False)
        #words, dist = flann.nn_index( descriptors.astype('float32') )
        
        for i,(c,ellipse) in enumerate(selected_mser):
            # Compute histogram
            hist = np.zeros(N, dtype='float')
            for ki,k in enumerate(keypoints):
                if is_in_ellipse(k.pt,ellipse):
                    hist[words[ki]] += 1

            # Normalize histogram
            norm = np.linalg.norm(hist)
            if norm > 0:
                hist /= norm

            cl = svc.predict(hist).flatten()
            
            if DEBUG:
                if cl == 1:
                    opacity = 0.4
                    img_color = cv2.cvtColor( img[z,:,:], cv2.cv.CV_GRAY2RGB )
                    for p in c:
                        img_color[p[1],p[0],:] = (
                            (1-opacity)*img_color[p[1],p[0],:]
                            + opacity * np.array([0,255,0])
                            )
                    cv2.imwrite(output_folder + "/"+str(z) + '_' +str(i) +"_region.png",img_color)

                img_color = cv2.cvtColor( img[z,:,:], cv2.cv.CV_GRAY2RGB )
                
                cv2.ellipse( img_color, (ellipse[0],
                                         (ellipse[1][0],ellipse[1][1]),
                                         ellipse[2]) , (0,0,255))
                for k_id,k in enumerate(keypoints):
                    if is_in_ellipse(k.pt,ellipse):
                        if cl == 1:
                            cv2.circle( img_color,
                                        (int(k.pt[0]),int(k.pt[1])),
                                        2,
                                        (0,255,0),
                                        -1)
                        else:
                            cv2.circle( img_color,
                                        (int(k.pt[0]),int(k.pt[1])),
                                        2,
                                        (0,0,255),
                                        -1)
                cv2.imwrite(output_folder + "/"+str(z) + '_' +str(i) +".png",img_color)

                cv2.ellipse( img_color_mser, (ellipse[0],
                                         (ellipse[1][0],ellipse[1][1]),
                                         ellipse[2]),
                             (0,255,0) if cl == 1 else (0,0,255) )

            if cl == 1:
                ellipse_center = [ellipse[0][0],ellipse[0][1],z]
                print c.shape
                if return_image_regions:
                    image_regions[-1].append((ellipse_center,c))
                else:
                    region = np.hstack( (c, [[z]]*c.shape[0]) )
                    detections[-1].append( (img.ImageToWorld(ellipse_center),
                                            img.ImageToWorld(region)) )

                
        if DEBUG:
            cv2.imwrite(output_folder + "/"+str(z) + "_color_mser.png",img_color_mser)

    if return_image_regions:
        return image_regions
    else:
        return detections
        
######################################################################

class BoxModel:
    def __init__(self,ofd,debug=False):
        self.debug = debug
        self.ofd = ofd
        
    def fit(self, data):
        # http://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
        mean = data[:,:3].mean(axis=0)
        # Do an SVD on the mean-centered data
        # u, s, v = np.linalg.svd(data - mean)
        # Now vv[0] contains the first principal component, i.e. the direction
        # vector of the 'best fit' line in the least squares sense.
        #return mean, v[0]
        return mean, np.array([0,0,1.0],dtype='float'), self.ofd
    
    def get_error( self, data, model):
        mean, u, ofd = model
        err_per_point = np.zeros(len(data),dtype='float')
        for i in range(len(data)):
            tmp = data[i,:3] - mean
            if np.linalg.norm(tmp) > ofd/2:
                err_per_point[i] = np.inf
            elif ( abs(data[i,3] - mean[0]) > ofd/2
                   or abs(data[i,4] - mean[0]) > ofd/2
                   or abs(data[i,5] - mean[1]) > ofd/2
                   or abs(data[i,6] - mean[1]) > ofd/2 ):
                err_per_point[i] = np.inf
            else:
                err_per_point[i] = np.linalg.norm(
                    tmp - np.dot(tmp,u)*u)
        
        return err_per_point
    
class SphereModel:
    def __init__(self,ofd,debug=False):
        self.debug = debug
        self.ofd = ofd
        
    def fit(self, data):
        mean = data[:,:3].mean(axis=0)
        return mean, self.ofd
    
    def get_error( self, data, model):
        mean, ofd = model
        err_per_point = np.zeros(len(data),dtype='float')
        for i in range(len(data)):
            distance = np.linalg.norm(data[i,:3] - mean)
            if distance > ofd/2:
                err_per_point[i] = np.inf
            elif ( data[i,3] < mean[0] - ofd/2
                   or data[i,4] > mean[0] + ofd/2
                   or data[i,5] < mean[1] - ofd/2
                   or data[i,6] > mean[1] + ofd/2
                   or data[i,7] < mean[2] - ofd/2
                   or data[i,8] > mean[2] + ofd/2 ):
                err_per_point[i] = np.inf
            else:
                err_per_point[i] = distance
        
        return err_per_point

def ransac_ellipses( detections,
                     ga,
                     model="sphere",
                     nb_iterations=100,
                     return_indices=False ):
    """
    RANSAC is performed in world coordinates
    """
    OFD = get_OFD(ga)

    centers = []
    for ellipse_center, region in detections:
        centers.append([ellipse_center[0], # x
                        ellipse_center[1], # y
                        ellipse_center[2], # z
                        region[:,0].min(),region[:,0].max(),
                        region[:,1].min(),region[:,1].max(),  
                        region[:,2].min(),region[:,2].max()])  

    centers = np.array(centers,dtype='float')
    if model == "sphere":
        model = SphereModel(OFD,debug=True)
    else:
        model = BoxModel(OFD,debug=True)
        
    # run RANSAC algorithm
    ransac_fit, ransac_data = ransac.ransac( centers,model,
                                             5, nb_iterations, 10.0, 10, # misc. parameters
                                             debug=True,return_all=True)

    if ransac_fit is None:
        print "RANSAC fiting failed"
        exit(1)

    inliers = ransac_data['inliers']

    if return_indices:
        return ransac_fit, inliers

    selected_regions = []
    for i,(ellipse_center, region) in enumerate(detections):
        if i in inliers:
            selected_regions.extend(region)

    return ransac_fit, selected_regions
