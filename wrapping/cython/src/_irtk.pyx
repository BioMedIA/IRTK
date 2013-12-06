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

from libcpp.vector cimport vector
from libcpp cimport bool

np.import_array()

ctypedef unsigned char uchar
ctypedef unsigned short ushort
ctypedef unsigned int uint

cdef extern from "irtk2cython.h":
     void Initialize()

def initialise_library():
    Initialize()

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
    void _transform_points( double* m,
                            double* pts,
                            size_t n )
    void _points_to_image( unsigned char* img,
                       double* pixelSize,
                       double* xAxis,
                       double* yAxis,
                       double* zAxis,
                       double* origin,
                       int* dim,
                       double* pts,
                       size_t n )
    void _orientation( double* pixelSize,
                       double* xAxis,
                       double* yAxis,
                       double* zAxis,
                       double* origin,
                       int* dim,
                       int &i, int &j, int &k)

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

    if dtype < 0:
        return None, None
        
    header = dict()
    header['pixelSize'] = pixelSize
    header['orientation'] = np.array( [xAxis, yAxis, zAxis],
                                      dtype="float64" )
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
    new_x = max(1,new_x)
    new_y = max(1,new_y)
    new_z = max(1,new_z)

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

def transform_points( np.ndarray[double, ndim=2,  mode="c"] m,
                      np.ndarray[double, ndim=2,  mode="c"] pts ):
    cdef size_t n = pts.shape[0]
    _transform_points( <double*> m.data,
                        <double*> pts.data,
                        n )
    return pts

def points_to_image( np.ndarray[double, ndim=2,  mode="c"] pts,
                     header ):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']

    cdef size_t n = pts.shape[0]

    cdef np.ndarray[unsigned char, ndim=4,  mode="c"] img = np.zeros( (dim[3],
                                                               dim[2],
                                                               dim[1],
                                                               dim[0]),
                                                              dtype='uint8' )

    _points_to_image( <unsigned char*> img.data,
                       <double*> pixelSize.data,
                       <double*> xAxis.data,
                       <double*> yAxis.data,
                       <double*> zAxis.data,
                       <double*> origin.data,
                       <int*> dim.data,
                       <double*> pts.data,
                       n )

    return img

def orientation( header ):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']

    orientation_code = { 1: "L2R", # Left to Right
                         2: "R2L", # Right to Left         
                         3: "P2A", # Posterior to Anterior
                         4: "A2P", # Anterior to Posterior
                         5: "I2S", # Inferior to Superior
                         6: "S2I " # Superior to Inferior
                         }
    
    cdef int i, j, k
    _orientation( <double*> pixelSize.data,
                   <double*> xAxis.data,
                   <double*> yAxis.data,
                   <double*> zAxis.data,
                   <double*> origin.data,
                   <int*> dim.data,
                   i, j, k )
    
    return ( orientation_code[i],
             orientation_code[j],
             orientation_code[k] )

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
                     uchar* img,
                     double* pixelSize,
                     double* xAxis,
                     double* yAxis,
                     double* zAxis,
                     double* origin,
                     int* dim )
    void _shrinkDisk( uchar* img,
                  int shape0,
                  int shape1,
                  double* center,
                  double radius,
                  int steps )

def voxellise( np.ndarray[double, ndim=2,  mode="c"] points,
               np.ndarray[int, ndim=2,  mode="c"] triangles,
               header ):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']

    cdef np.ndarray[uchar, ndim=4,  mode="c"] img = np.zeros( (dim[3],
                                                               dim[2],
                                                               dim[1],
                                                               dim[0]),
                                                              dtype='uint8' )

    _voxellise( <double*> points.data,
                 points.shape[0],
                 <int*> triangles.data,
                 triangles.shape[0],
                 <uchar*> img.data,
                 <double*> pixelSize.data,
                 <double*> xAxis.data,
                 <double*> yAxis.data,
                 <double*> zAxis.data,
                 <double*> origin.data,
                 <int*> dim.data )

    return img

def shrinkDisk( np.ndarray[uchar, ndim=2,  mode="c"] img,
               np.ndarray[double, ndim=1,  mode="c"] center,
               double radius,
                int steps ):
    _shrinkDisk( <uchar*> img.data,
                  img.shape[0],
                  img.shape[1],
                  <double*> center.data,
                  radius,
                  steps )
    return img


########## Functions on list of images ##########

cdef extern from "irtk2cython.h":
    void _write_list( float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,
                 int n )

def write_list( np.ndarray[float, ndim=1,  mode="c"] img,
                np.ndarray[double, ndim=1,  mode="c"] pixelSize,
                np.ndarray[double, ndim=1,  mode="c"] xAxis,
                np.ndarray[double, ndim=1,  mode="c"] yAxis,
                np.ndarray[double, ndim=1,  mode="c"] zAxis,
                np.ndarray[double, ndim=1,  mode="c"] origin,
                np.ndarray[int, ndim=1,  mode="c"] dim,
                int n ):
    _write_list( <float*> img.data,
                 <double*> pixelSize.data,
                 <double*> xAxis.data,
                 <double*> yAxis.data,
                 <double*> zAxis.data,
                 <double*> origin.data,
                 <int*> dim.data,
                 n)


########## Reconstruction ##########
cdef extern from "reconstruction.h":
    void _reconstruct(
                 # input stacks or slices
                 float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,

                 # number of stacks
                 int n,

                 # stack ids: which stack each slice
                 # comes from
                 int* _stack_ids,

                 # number of reconstruction iterations to run
                 int iterations,

                 # initial transformations
                 double* tx,
                 double* ty,
                 double* tz,
                 double* rx,
                 double* ry,
                 double* rz,

                 # slice thickness
                 double* _thickness,

                 # mask (header same as template)
                 float* mask_img,
                 
                 # output: reconstructed image
                 float* reconstructed_img,
                 double* reconstructed_pixelSize,
                 double* reconstructed_xAxis,
                 double* reconstructed_yAxis,
                 double* reconstructed_zAxis,
                 double* reconstructed_origin,
                 int* reconstructed_dim )

def reconstruct( np.ndarray[float, ndim=1,  mode="c"] img,
                 np.ndarray[double, ndim=1,  mode="c"] pixelSize,
                 np.ndarray[double, ndim=1,  mode="c"] xAxis,
                 np.ndarray[double, ndim=1,  mode="c"] yAxis,
                 np.ndarray[double, ndim=1,  mode="c"] zAxis,
                 np.ndarray[double, ndim=1,  mode="c"] origin,
                 np.ndarray[int, ndim=1,  mode="c"] dim,
                 int n,
                 np.ndarray[int, ndim=1,  mode="c"] stack_ids,
                 int iterations,
                 np.ndarray[double, ndim=1,  mode="c"] tx,
                 np.ndarray[double, ndim=1,  mode="c"] ty,
                 np.ndarray[double, ndim=1,  mode="c"] tz,
                 np.ndarray[double, ndim=1,  mode="c"] rx,
                 np.ndarray[double, ndim=1,  mode="c"] ry,
                 np.ndarray[double, ndim=1,  mode="c"] rz,
                 np.ndarray[double, ndim=1,  mode="c"] thickness,
                 np.ndarray[float, ndim=4,  mode="c"] mask,
                 template_header ):

    cdef np.ndarray[double, ndim=1,  mode="c"] template_pixelSize = template_header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] template_xAxis = template_header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] template_yAxis = template_header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] template_zAxis = template_header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] template_origin = template_header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] template_dim =  template_header['dim']

    cdef np.ndarray[float, ndim=4,  mode="c"] reconstructed = np.zeros( (template_dim[3],
                                                                         template_dim[2],
                                                                         template_dim[1],
                                                                         template_dim[0]),
                                                                        dtype='float32' )

    _reconstruct(
        # input stacks or slices
        <float*> img.data,
         <double*> pixelSize.data,
         <double*> xAxis.data,
         <double*> yAxis.data,
         <double*> zAxis.data,
         <double*> origin.data,
         <int*> dim.data,

         # number of stacks
         n,

         # stack ids: which stack each slice
         # comes from
         <int*>stack_ids.data,

         # number of reconstruction iterations to run
         iterations,

         # initial transformations
         <double*> tx.data,
         <double*> ty.data,
         <double*> tz.data,
         <double*> rx.data,
         <double*> ry.data,
         <double*> rz.data,

         # slice thickness
         <double*> thickness.data,

         # mask (header same as template)
         <float*> mask.data,
         
         # output: reconstructed image
         <float*> reconstructed.data,
         <double*> template_pixelSize.data,
         <double*> template_xAxis.data,
         <double*> template_yAxis.data,
         <double*> template_zAxis.data,
         <double*> template_origin.data,
         <int*> template_dim.data )

    return reconstructed

########## Drawing ##########

cdef extern from "drawing.h":
    void _drawSphere( unsigned char* img,
                  int shape0,
                  int shape1,
                  int shape2,
                  int x0,
                  int y0,
                  int z0,
                  int rx,
                  int ry,
                  int rz )

def drawSphere( np.ndarray[unsigned char, ndim=3,  mode="c"] img,
                int x0, int y0, int z0,
                int rx, int ry, int rz ):
    cdef int shape0 = img.shape[0]
    cdef int shape1 = img.shape[1]
    cdef int shape2 = img.shape[2]
    _drawSphere( <unsigned char*> img.data,
                  shape0, shape1, shape2,
                  x0, y0, z0,
                  rx, ry, rz )
    print "MAX:", np.max(img)
    return img

#### Points ####

cdef extern from "irtk2cython.h":
    void _read_points( char *filename,
                       vector[ vector[double] ] &points )
    void _write_points( char *filename,
                        vector[ vector[double] ] &points )

def read_points( bytes py_string ):
    cdef char* c_string = py_string
    cdef vector[vector[double]] v
    _read_points( c_string, v )
    return v

cdef extern from "registration.h":
    double _registration_rigid_points( double* source_points,
                                       double* target_points,
                                       int n,
                                       double &tx,
                                       double &ty,
                                       double &tz,
                                       double &rx,
                                       double &ry,
                                       double &rz )

def registration_rigid_points( np.ndarray[double, ndim=2,  mode="c"] source,
                               np.ndarray[double, ndim=2,  mode="c"] target ):
    cdef double tx
    cdef double ty
    cdef double tz
    cdef double rx
    cdef double ry
    cdef double rz
    cdef double RMS = _registration_rigid_points( <double*> source.data,
                                                   <double*> target.data,
                                                   source.shape[0],
                                                   tx, ty, tz,
                                                   rx, ry,rz )
    return ( ( tx, ty, tz, rx, ry, rz ), RMS )

### CRF ###

ctypedef double pixel_t

ctypedef short LabelID
ctypedef long long EnergyType

cdef extern from "crf.h":
    void _crf( pixel_t* img,
               double* pixelSize,
               double* xAxis,
               double* yAxis,
               double* zAxis,
               double* origin,
               int* dim,
               LabelID* labels,
               double* proba,
               double l,
               double sigma,
               double sigmaZ )
    
def crf( np.ndarray[pixel_t, ndim=4, mode="c"] img,
           header,
           np.ndarray[LabelID, ndim=4, mode="c"] labels,
           np.ndarray[double, ndim=4, mode="c"] proba,
           double l=1.0,
           double sigma=0.0,
           double sigmaZ=0.0 ):
    cdef np.ndarray[double, ndim=1,  mode="c"] pixelSize = header['pixelSize']
    cdef np.ndarray[double, ndim=1,  mode="c"] xAxis = header['orientation'][0]
    cdef np.ndarray[double, ndim=1,  mode="c"] yAxis = header['orientation'][1]
    cdef np.ndarray[double, ndim=1,  mode="c"] zAxis = header['orientation'][2]
    cdef np.ndarray[double, ndim=1,  mode="c"] origin = header['origin']
    cdef np.ndarray[int, ndim=1,  mode="c"] dim =  header['dim']    
    print "starting crf..."                                                                 
    _crf( <pixel_t*> img.data,
           <double*> pixelSize.data,
           <double*> xAxis.data,
           <double*> yAxis.data,
           <double*> zAxis.data,
           <double*> origin.data,
           <int*> dim.data,
           <LabelID*> labels.data,
           <double*> proba.data,
           l, sigma, sigmaZ )
    
    return labels

