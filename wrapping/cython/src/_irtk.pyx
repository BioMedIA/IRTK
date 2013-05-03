#clib segmentation++
#clib registration2++
#clib registration++
#clib transformation++
#clib contrib++
#clib image++
#clib geometry++
#clib common++
#clib recipes
#clib znz
#clib niftiio

import numpy as np
cimport numpy as np

np.import_array()

ctypedef unsigned char uchar
ctypedef unsigned short ushort
ctypedef unsigned int uint

########## Image function ##########

cdef extern from "irtk2cython.h":
    int _get_header( char* filename,
                      double* pixelSize,
                      double* xAxis,
                      double* yAxis,
                      double* zAxis,
                      double* origin,
                      int* dim )
    void _resample( float* img_in,
                double* pixelSize,
                double* xAxis,
                double* yAxis,
                double* zAxis,
                double* origin,
                int* dim,
                double* new_pixelSize,
                float* img_out,
                int interpolation_method,
                float gaussian_parameter )
    
def get_header( bytes py_string ):
    cdef char* c_string = py_string

    types = {
        1: "int8",
        2: "uint8",
        3: "int16",
        4: "uint16",
        5: "int32",
        6: "uint32",
        7: "float32",
        8: "float64"
        }
    
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = np.zeros( 4, dtype="float64" )
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = np.zeros( 3, dtype="float64" )
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = np.zeros( 3, dtype="float64" )
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = np.zeros( 3, dtype="float64" )
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = np.zeros( 4, dtype="float64" )
    cdef np.ndarray[int, ndim=1,  mode="c"] dim = np.zeros( 4, dtype="int32" )

    cdef int dtype = _get_header( c_string,
                                  <double*> pixelSize.data,
                                  <double*> xAxis.data,
                                  <double*> yAxis.data,
                                  <double*> zAxis.data,
                                  <double*> origin.data,
                                  <int*> dim.data )

    header = dict()
    header['pixelSize'] = pixelSize
    header['orientation'] = (xAxis, yAxis, zAxis)
    header['origin'] = origin
    header['dim'] = dim

    return header, types[dtype]

def resample( np.ndarray[float, ndim=4,  mode="c"] img_in,
              header,
              np.ndarray[double, ndim=1,  mode="c"] new_pixelSize,
              interpolation_method='linear',
              float gaussian_parameter=1.0):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']

    # Determine the new dimensions of the image
    cdef int new_x = int(float(dim[0]) * pixelSize[0] / new_pixelSize[0])
    cdef int new_y = int(float(dim[1]) * pixelSize[1] / new_pixelSize[1])
    cdef int new_z = int(float(dim[2]) * pixelSize[2] / new_pixelSize[2])

    cdef np.ndarray[float, ndim=4,  mode="c"] img_out = np.zeros( (dim[3], # T
                                                                   new_z,
                                                                   new_y,
                                                                   new_x),
                                                                  dtype='float32' )

    interpolation = { 'nearest' : 0,
                      'linear' : 1,
                      'bspline' : 2,
                      'cspline' : 3,
                      'sinc' : 4,
                      'shape' : 5,
                      'gaussian' : 6
                      }

    cdef int inter = interpolation[interpolation_method]

    
    _resample( <float*> img_in.data,
                <double*> pixelSize.data,
                <double*> xAxis.data,
                <double*> yAxis.data,
                <double*> zAxis.data,
                <double*> origin.data,
                <int*> dim.data,
                <double*> new_pixelSize.data,
                <float*> img_out.data,
                inter,
                gaussian_parameter )

    return img_out, header
    
########## Registration ##########

cdef extern from "registration.h":
    void _read_rigid( char* filename,
                      double &tx,
                      double &ty,
                      double &tz,
                      double &rx,
                      double &ry,
                      double &rz )
    void _write_rigid( char* filename,
                       double tx,
                       double ty,
                       double tz,
                       double rx,
                       double ry,
                       double rz )
    void _transform_rigid( double tx,
                           double ty,
                           double tz,
                           double rx,
                           double ry,
                           double rz,
                           float* source_img,
                           double* source_pixelSize,
                           double* source_xAxis,
                           double* source_yAxis,
                           double* source_zAxis,
                           double* source_origin,
                           int* source_dim,
                           float* target_img,
                           double* target_pixelSize,
                           double* target_xAxis,
                           double* target_yAxis,
                           double* target_zAxis,
                           double* target_origin,
                           int* target_dim,
                           int interpolation_method,
                           float gaussian_parameter )
    void _registration_rigid( short* source_img,
                              double* source_pixelSize,
                              double* source_xAxis,
                              double* source_yAxis,
                              double* source_zAxis,
                              double* source_origin,
                              int* source_dim,
                              short* target_img,
                              double* target_pixelSize,
                              double* target_xAxis,
                              double* target_yAxis,
                              double* target_zAxis,
                              double* target_origin,
                              int* target_dim,
                              double &tx,
                              double &ty,
                              double &tz,
                              double &rx,
                              double &ry,
                              double &rz )

def read_rigid( bytes py_string ):
    cdef char* c_string = py_string
    cdef double tx, ty, tz, rx, ry, rz
    _read_rigid( c_string,
                 tx, ty, tz,
                 rx, ry, rz )
    return (tx, ty, tz, rx, ry, rz )

def write_rigid( bytes py_string,
                 double tx,
                 double ty,
                 double tz,
                 double rx,
                 double ry,
                 double rz ):
    cdef char* c_string = py_string
    _write_rigid( c_string,
                  tx, ty, tz,
                  rx, ry, rz )
    return

def transform_rigid( double tx,
                     double ty,
                     double tz,
                     double rx,
                     double ry,
                     double rz,
                     np.ndarray[float, ndim=4,  mode="c"] source_img,
                     source_header,
                     target_header,
                     interpolation_method,
                     float gaussian_parameter ):
    cdef np.ndarray[double, ndim=1,  mode="c"] source_pixelSize = source_header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] source_xAxis = source_header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] source_yAxis = source_header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] source_zAxis = source_header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] source_origin = source_header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] source_dim =  source_header['dim']

    cdef np.ndarray[double, ndim=1,  mode="c"] target_pixelSize = target_header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] target_xAxis = target_header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] target_yAxis = target_header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] target_zAxis = target_header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] target_origin = target_header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] target_dim =  target_header['dim']

    cdef np.ndarray[float, ndim=4,  mode="c"] target_img = np.zeros( (target_dim[3],
                                                                      target_dim[2],
                                                                      target_dim[1],
                                                                      target_dim[0]),
                                                                     dtype='float32' )

    interpolation = { 'nearest' : 0,
                      'linear' : 1,
                      'bspline' : 2,
                      'cspline' : 3,
                      'sinc' : 4,
                      'shape' : 5,
                      'gaussian' : 6
                      }

    cdef int inter = interpolation[interpolation_method]
    
    _transform_rigid( tx, ty, tz,
                      rx, ry,rz,
                      <float*> source_img.data,
                      <double*> source_pixelSize.data,
                      <double*> source_xAxis.data,
                      <double*> source_yAxis.data,
                      <double*> source_zAxis.data,
                      <double*> source_origin.data,
                      <int*> source_dim.data,
                      <float*> target_img.data,
                      <double*> target_pixelSize.data,
                      <double*> target_xAxis.data,
                      <double*> target_yAxis.data,
                      <double*> target_zAxis.data,
                      <double*> target_origin.data,
                      <int*> target_dim.data,
                      inter,
                      gaussian_parameter )

    return target_img

def registration_rigid( np.ndarray[short, ndim=4,  mode="c"] source_img,
                        source_header,
                        np.ndarray[short, ndim=4,  mode="c"] target_img,
                        target_header,
                        double tx,
                        double ty,
                        double tz,
                        double rx,
                        double ry,
                        double rz ):
    cdef np.ndarray[double, ndim=1,  mode="c"] source_pixelSize = source_header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] source_xAxis = source_header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] source_yAxis = source_header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] source_zAxis = source_header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] source_origin = source_header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] source_dim =  source_header['dim']

    cdef np.ndarray[double, ndim=1,  mode="c"] target_pixelSize = target_header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] target_xAxis = target_header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] target_yAxis = target_header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] target_zAxis = target_header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] target_origin = target_header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] target_dim =  target_header['dim']

    _registration_rigid( <short*> source_img.data,
                      <double*> source_pixelSize.data,
                      <double*> source_xAxis.data,
                       <double*> source_yAxis.data,
                       <double*> source_zAxis.data,
                       <double*> source_origin.data,
                       <int*> source_dim.data,
                       <short*> target_img.data,
                       <double*> target_pixelSize.data,
                       <double*> target_xAxis.data,
                       <double*> target_yAxis.data,
                       <double*> target_zAxis.data,
                       <double*> target_origin.data,
                       <int*> target_dim.data,
                       tx, ty, tz,
                       rx, ry,rz )

    return ( tx, ty, tz, rx, ry, rz )

########## Voxellise ##########

cdef extern from "voxellise.h":
    void _voxellise( double* points,
                     int npoints,
                     int* triangles,
                     int ntriangles,
                     short* img,
                     double* pixelSize,
                     double* xAxis,
                     double* yAxis,
                     double* zAxis,
                     double* origin,
                     int* dim )

def voxellise( np.ndarray[double, ndim=2,  mode="c"] points,
               np.ndarray[int, ndim=2,  mode="c"] triangles,
               header ):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']

    cdef np.ndarray[short, ndim=4,  mode="c"] img = np.zeros( (dim[3],
                                                               dim[2],
                                                               dim[1],
                                                               dim[0]),
                                                              dtype='int16' )

    _voxellise( <double*> points.data,
                 points.shape[0],
                 <int*> triangles.data,
                 triangles.shape[0],
                 <short*> img.data,
                 <double*> pixelSize.data,
                 <double*> xAxis.data,
                 <double*> yAxis.data,
                 <double*> zAxis.data,
                 <double*> origin.data,
                 <int*> dim.data )

    return img
