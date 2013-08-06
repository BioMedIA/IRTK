#!/usr/bin/python

import cv2
import sys
import csv
import numpy as np
from math import cos,sin
import math
import SimpleITK as sitk
import argparse

from glob import glob
from sklearn import cluster
from sklearn import neighbors
from sklearn import svm
from sklearn.externals import joblib
from joblib import Parallel, delayed

from scipy.stats.mstats import mquantiles

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

def process_file( raw_file, ga, coordinates, size, classifier, N, NEW_SAMPLING, DEBUG=False ):
    X = []
    Y = []

    ofd_model = np.array([ -4.97315445e-03,
                            3.19846853e-01,
                            -2.60839214e+00,
                            2.62679565e+01])

    OFD = 0
    for k in range(4):
        OFD += ofd_model[3-k]*ga**k

    OFD /= NEW_SAMPLING

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

    coordinates = np.array(map(float,coordinates.split(',')),dtype='float')
    size = np.array(map(float,size.split(',')),dtype='float')

    sitk_img = sitk.ReadImage( raw_file )
    raw_spacing = sitk_img.GetSpacing()

    ## Resample
    resample = sitk.ResampleImageFilter()
    resample.SetOutputDirection(sitk_img.GetDirection())
    resample.SetOutputOrigin(sitk_img.GetOrigin())
    resample.SetOutputSpacing([NEW_SAMPLING,NEW_SAMPLING,raw_spacing[2]])
    resample.SetSize([int(sitk_img.GetSize()[0]*raw_spacing[0]/NEW_SAMPLING),
                      int(sitk_img.GetSize()[1]*raw_spacing[1]/NEW_SAMPLING),
                      sitk_img.GetSize()[2]])
    sitk_img = resample.Execute(sitk_img)

    # Adjust coordinates and size
    box = coordinates * np.array([1.0,raw_spacing[1]/NEW_SAMPLING,raw_spacing[0]/NEW_SAMPLING],dtype='float')
    box_size = size * np.array([1.0,raw_spacing[1]/NEW_SAMPLING,raw_spacing[0]/NEW_SAMPLING],dtype='float')
    z0,y0,x0 = box.astype('int')
    d0,h0,w0 = box_size.astype('int')

    brain_center = (x0 + w0/2, y0 + h0/2)
    
    data = sitk.GetArrayFromImage( sitk_img ).astype("float")

    ## Contrast-stretch with saturation
    q = mquantiles(data.flatten(),[0.01,0.99])
    data[data<q[0]] = q[0]
    data[data>q[1]] = q[1]
    data -= data.min()
    data /= data.max()
    data *= 255
    data = data.astype('uint8')

    for z in range(data.shape[0]):
        contours = mser.detect(data[z,:,:])
        keypoints = sift.detect(data[z,:,:])
        if keypoints is None or len(keypoints) == 0:
            continue
        (keypoints, descriptors) = siftExtractor.compute(data[z,:,:],keypoints)     

        for i,c in enumerate(contours):
            hist = np.zeros(N, dtype='float')
            ellipse = cv2.fitEllipse(np.array(map(lambda x:[x],
                                                  c),dtype='int32'))
            
            # filter by size
            if ( ellipse[1][0] > OFD
                 or ellipse[1][1] > OFD
                 or ellipse[1][0] < 0.5*OFD
                 or ellipse[1][1] < 0.5*OFD ) :
                continue
            
            # filter by eccentricity
            # if math.sqrt(1-(np.min(ellipse[1])/np.max(ellipse[1]))**2) > 0.75:
            #     continue

            distance = math.sqrt((ellipse[0][0]-brain_center[0])**2
                                 +(ellipse[0][1]-brain_center[1])**2)

            if max(w0,h0)/2 >= distance >= min(w0,h0)/8:
                continue            

            for k,d in zip(keypoints,descriptors):
                if is_in_ellipse(k.pt,ellipse):
                    c = classifier.kneighbors(d, return_distance=False)
                    hist[c] += 1

            # Normalize histogram
            norm = np.linalg.norm(hist)
            if norm > 0:
                hist /= norm

            if distance > max(w0,h0)/4:
                if DEBUG: print 0
                X.append(hist)
                Y.append(0)
            else:
                if distance < min(w0,h0)/8 and z0 + d0/8 <= z <= z0+7*d0/8:
                    if DEBUG: print 1
                    X.append(hist)
                    Y.append(1)
                else:
                    continue
                   
            if DEBUG:
                img_color = cv2.cvtColor( data[z,:,:], cv2.cv.CV_GRAY2RGB )
                cv2.ellipse( img_color, (ellipse[0],
                                         (ellipse[1][0],ellipse[1][1]),
                                         ellipse[2]) , (0,0,255))
                for k_id,k in enumerate(keypoints):
                    if is_in_ellipse(k.pt,ellipse):
                        if Y[-1] == 1:
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
                cv2.imwrite("/tmp/"+str(z) + '_' +str(i) +'_'+str(k_id)+".png",img_color)
                # cv2.imshow("show",img_color)
                # cv2.waitKey(0)
    return X,Y
        
######################################################################

parser = argparse.ArgumentParser(
    description='Learn MSER classifier using SIFT BOW.' )
parser.add_argument( '--training_patients' )
parser.add_argument( '--original_folder' )
parser.add_argument( '--ga_file' )
parser.add_argument( '--clean_brainboxes' )
parser.add_argument( '--new_sampling', type=float )
parser.add_argument( '--vocabulary' )
parser.add_argument( '--output' )
parser.add_argument( '--debug', action="store_true", default=False )

args = parser.parse_args()
    
vocabulary = open(args.vocabulary, 'rb')
voca = np.load(vocabulary)
classifier = neighbors.NearestNeighbors(1)
N = voca.shape[0] 
classifier.fit(voca)

f = open( args.training_patients, "r" )
patients = []
for p in f:
    patients.append(p.rstrip())
f.close()

reader = csv.reader( open( args.ga_file, "rb"), delimiter=" " )
all_ga = {}
for patient_id, ga in reader:
    all_ga[patient_id] = float(ga)
    
reader = csv.reader( open( args.clean_brainboxes, "r" ),
                     delimiter='\t' )

training_patients = []
for patient_id, raw_file, cl, coordinates, size in reader:
    if patient_id not in patients:
        print "Skipping testing patient: " + patient_id
        continue
    training_patients.append((args.original_folder + '/' + raw_file,all_ga[patient_id],coordinates, size))

XY = Parallel(n_jobs=-1)(delayed(process_file)(raw_file,ga,coordinates, size, classifier,N,args.new_sampling,args.debug)
                         for raw_file,ga,coordinates, size in training_patients )

print len(XY)

X = []
Y = []
for x,y in XY:
    X.extend(x)
    Y.extend(y)

print "RATIO = ", np.sum(Y), len(Y)
    
X = np.array(X,dtype='float')
Y = np.array(Y,dtype='float')

# svc = svm.SVC()

svc = svm.LinearSVC(dual=False,)

svc.fit(X, Y)
print svc.score(X, Y)

joblib.dump(svc, args.output)
