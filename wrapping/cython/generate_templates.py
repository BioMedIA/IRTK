#!/usr/bin/python

import sys

output_folder = sys.argv[1]

#define IRTK_VOXEL_UNKNOWN  				0
#define IRTK_VOXEL_CHAR     				1
#define IRTK_VOXEL_UNSIGNED_CHAR    2
#define IRTK_VOXEL_SHORT    				3
#define IRTK_VOXEL_UNSIGNED_SHORT   4
#define IRTK_VOXEL_INT      				5
#define IRTK_VOXEL_UNSIGNED_INT     6
#define IRTK_VOXEL_FLOAT    				7
#define IRTK_VOXEL_DOUBLE   				8
#define IRTK_VOXEL_RGB      				9

types = {
    "char" : "int8",
    "uchar" : "uint8",
    "short" : "int16",
    "ushort" : "uint16",
    "int" : "int32",
    "uint" : "uint32",
    "float" : "float32",
    "double" : "float64"
    }

header_cpp = {
"""void _imread_%s( char* filename, %s* img );""" : 2,
"""
void _imwrite_%s( char* filename,
                  %s* img,
                  double* pixelSize,
                  double* xAxis,
                  double* yAxis,
                  double* zAxis,
                  double* origin,
                  int* dim );
""" : 2,
}

templates_cpp = {
"""
void _imread_%s( char* filename,
                  %s* img ) {
    irtkGenericImage<%s> image(filename);
    int n   = image.GetNumberOfVoxels();
    %s* ptr = image.GetPointerToVoxels();
    for ( int i = 0; i < n; i++)
        img[i] = ptr[i];
}
""" : 4,

"""
void _imwrite_%s( char* filename,
                  %s* img,
                  double* pixelSize,
                  double* xAxis,
                  double* yAxis,
                  double* zAxis,
                  double* origin,
                  int* dim ) {
    irtkGenericImage<%s> irtk_image;
    py2irtk<%s>( irtk_image,
                 img,
                 pixelSize,
                 xAxis,
                 yAxis,
                 zAxis,
                 origin,
                 dim );
    irtk_image.Write( filename );
}
""" : 4,

}

header_py = {
"""void _imread_%s( char* filename, %s* img )""" : 2,
"""void _imwrite_%s( char* filename,
                  %s* img,
                  double* pixelSize,
                  double* xAxis,
                  double* yAxis,
                  double* zAxis,
                  double* origin,
                  int* dim )
""" : 2,
}

templates_py = {
"""
def imread_%s( bytes py_string, header ):
    cdef char* c_string = py_string
    cdef np.ndarray[%s, ndim=4,  mode="c"] img = np.empty( (header['dim'][3],
                                                            header['dim'][2],
                                                            header['dim'][1],
                                                            header['dim'][0]), dtype="%s" )
    _imread_%s( c_string, <%s*>img.data )
    return img
""" : (0,0,1,0,0),

"""
def imwrite_%s( bytes py_string, np.ndarray[%s, ndim=4,  mode="c"] img, header ):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']

    cdef char* c_string = py_string
    _imwrite_%s( c_string,
                <%s*>img.data,
                <double*> pixelSize.data,
                <double*> xAxis.data,
                <double*> yAxis.data,
                <double*> zAxis.data,
                <double*> origin.data,
                <int*> dim.data )
""" : (0,0,0,0),
}

f = open( output_folder + "/templates.h",'wb')
f.write("""#ifndef TEMPLATES_H
#define TEMPLATES_H

#include "irtk2cython.h"\n\n""")

for t,n in header_cpp.iteritems():
    for dtype in types:
        f.write( t % ((dtype,)*n) )
        f.write("\n")
    f.write("\n")
    
f.write("""#endif""")
f.close()

f = open( output_folder + "/templates.cc",'wb')
f.write("""#include "templates.h"\n""")

for t,n in templates_cpp.iteritems():
    for dtype in types:
        f.write( t % ((dtype,)*n) )

f.close()

f = open( output_folder + "/templates.pyx",'wb')
f.write(
 """########## Automatically generated templated code ##########

cdef extern from "templates.h":\n""")

for t,n in header_py.iteritems():
    for dtype in types:
        f.write('    ')
        f.write( t % ((dtype,)*n) )
        f.write("\n")
f.write("\n")

for t,n in templates_py.iteritems():
    for s in types.iteritems():
        n2 = map( lambda x: s[x], n )
        f.write( t % tuple(n2) )
        
f.close()
