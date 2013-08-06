# http://www.tp.umu.se/~nylen/fnm/pylect/advanced/image_processing/index.html
# http://step.polymtl.ca/~rv101/energy/

import numpy as np
import pycrf
import itertools
import scipy
import cv2

def index( i, j, img ):
    return i + img.shape[0]*j


# Create Noisy Lena
lena = cv2.imread("lena.png",0)
lena = lena[230:310, 210:350]
#noisy = lena
noisy = lena + 4.5*lena.std()*np.random.random(lena.shape)
noisy[noisy<noisy.mean()] = 0
noisy[noisy>0] = 1
cv2.imwrite("noisy_lena.png",noisy*255)

num_pixels = lena.shape[0]*lena.shape[1]
num_labels = 2
crf = pycrf.CRF( num_pixels, num_labels )

crf.setLabels(noisy.flatten().astype('int32'))

datacost = np.zeros(num_pixels*num_labels,dtype='float32')
for p in xrange(num_pixels):
    for l in xrange(num_labels):
        datacost[p*num_labels+l] = (noisy.flat[p]-l)**2
crf.setDataCost(datacost)
print len(datacost),datacost

labelcost = np.zeros(num_labels*num_labels,dtype='float32')
for l1 in xrange(num_labels):
    for l2 in xrange(num_labels):
        labelcost[l1+num_labels*l2] = int(l1 != l2)*1100
crf.setSmoothCost(labelcost) 
print labelcost

r = [-1,0, 1]
neighbourhood = []
for a,b in itertools.product(r,r):
    if a+b != 0:
    #if abs(a)+abs(b) == 1:
        neighbourhood.append((a,b))

print neighbourhood
        
for i in xrange(lena.shape[0]):
    for j in xrange(lena.shape[1]):
        for a,b in neighbourhood:
            if ( 0 <= i+a < lena.shape[0]
                 and 0 <= j+b < lena.shape[1]
                 and index(i,j,lena) < index(i+a,j+b,lena) ):
                crf.setNeighbors( index(i,j,lena),
                                  index(i+a,j+b,lena),
                                  1)

print crf.compute_energy()
# for l in range(100,2000,500):
#     crf.setSmoothCost(labelcost*l) 
crf.expansion()
print crf.compute_energy()

denoised = np.reshape( crf.getLabels()*255, lena.shape)
cv2.imwrite("denoised_lena.png",denoised)
