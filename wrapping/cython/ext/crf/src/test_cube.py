# http://www.tp.umu.se/~nylen/fnm/pylect/advanced/image_processing/index.html
# http://step.polymtl.ca/~rv101/energy/

import numpy as np
import pycrf
import itertools
import scipy
import cv2

def index( i, j, img ):
    return i + img.shape[0]*j

img = cv2.imread("img/restoration-observed.png",0)
noisy = img

num_pixels = img.shape[0]*img.shape[1]
num_labels = 256
crf = pycrf.CRF( num_pixels, num_labels )

crf.setLabels(noisy.flatten().astype('int32'))

datacost = np.zeros(num_pixels*num_labels,dtype='float32')
# for p in xrange(num_pixels):
#     for l in xrange(num_labels):
#         datacost[p*num_labels+l] = (noisy.flat[p]-l)**2
crf.setDataCost(datacost)
print len(datacost),datacost

labelcost = np.zeros(num_labels*num_labels,dtype='float32')
for l1 in xrange(num_labels):
    for l2 in xrange(num_labels):
        labelcost[l1+num_labels*l2] = int(l1 != l2)#*2331
        #labelcost[l1+num_labels*l2] = min(2,abs(l1 - l2))*2331
crf.setSmoothCost(labelcost) 
print labelcost

r = [-1,0, 1]
neighbourhood = []
for a,b in itertools.product(r,r):
    if abs(a)+abs(b) == 1:
        neighbourhood.append((a,b))

print neighbourhood
        
for i in xrange(img.shape[0]):
    for j in xrange(img.shape[1]):
        for a,b in neighbourhood:
            if ( 0 <= i+a < img.shape[0]
                 and 0 <= j+b < img.shape[1]
                 and index(i,j,img) < index(i+a,j+b,img) ):
                crf.setNeighbors( index(i,j,img),
                                  index(i+a,j+b,img),
                                  1)

print crf.compute_energy()
crf.expansion()
print crf.compute_energy()

denoised = np.reshape( crf.getLabels(), img.shape)
cv2.imwrite("img/denoised_cube.png",denoised)
