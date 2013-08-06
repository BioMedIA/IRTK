#!/usr/bin/python

import sys
import numpy as np
import irtk

full_file = sys.argv[1]
cropped_file = sys.argv[2]

full_img = irtk.imread( full_file, dtype='float32' )
cropped_img = irtk.imread( cropped_file, dtype='float32' )

(z,y,x), score = irtk.match_template( full_img, cropped_img, pad_input=False )

print score

print ' '.join( map(str, [ full_img.shape[0],
                           full_img.shape[1],
                           full_img.shape[2],
                           z,
                           y,
                           x ] ) )

