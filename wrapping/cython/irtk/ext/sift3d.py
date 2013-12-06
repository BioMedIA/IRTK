#!/usr/bin/python

import csv
import numpy as np
import sys
import os
import subprocess
import random
import irtk

def sift3d( img, file_id=None, tmp_dir='tmp' ):

    if file_id is None:
        file_id = str( random.randint(1, 1000) )
        
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    tmp_file = tmp_dir + '/' + str(os.getpid()) + "_" + file_id + ".nii"
    tmp_sift = tmp_dir + '/' + str(os.getpid()) + "_" + file_id + ".txt"
    irtk.imwrite( tmp_file, img )

    proc = subprocess.Popen([os.path.dirname(irtk.__file__) + "/ext/featExtract.ubu",
                         tmp_file,
                         tmp_sift], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()

    print out
    print err

    os.remove( tmp_file )

    f = open(tmp_sift, "r")

    # skipping header
    line1 = f.next()
    line2 = f.next()
    line3 = f.next()

    titles = ["x", "y", "z", "scale", # Scale-space location
              "o11", "o12", "o13", "o21", "o22", "o23", "o31", "o32", "o32", # orientation
              "e1", "e2", "e3", # 2nd moment eigenvalues
              "i1" # info flag
              ]

    # descriptor
    for i in range(64):
        titles.append("d"+str(i))

    reader = csv.DictReader(f, titles, delimiter='\t',quoting=csv.QUOTE_NONNUMERIC)

    features = []
    for row in reader:
        descr = []
        for i in range(64):
            descr.append(float(row["d"+str(i)]))
        descr = np.array( descr, dtype='float' )
        features.append( ( ( int(row['x']),
                             int(row['y']),
                             int(row['z']) ),
                           descr )
                         )

    return features
