#!/usr/bin/python

import sys

import numpy as np
import irtk

# img1 = irtk.imread("reconstruction/t2.nii", dtype='float32')
# img2 = irtk.imread("reconstruction/t2seg.nii", dtype='float32')
# irtk.imwrite( "reconstruction/segfixed.nii", img2.transform(target=img1) )

shape = map( int, sys.argv[1:4] )
z = int(sys.argv[4])
y = int(sys.argv[5])
x = int(sys.argv[6])

target = irtk.imread( sys.argv[7] )

img = irtk.imread( sys.argv[8], dtype='float32' )

new_data = np.zeros( shape, dtype='int32' )
new_data[z:z+img.shape[0],
         y:y+img.shape[1],
         x:x+img.shape[2]] = img

new_img = irtk.Image( new_data, target.get_header() )
irtk.imwrite( sys.argv[9], new_img ) 
