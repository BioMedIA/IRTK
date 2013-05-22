#ifndef IRTK2CYTHON_H
#define IRTK2CYTHON_H

#include <string>
#include <sstream>
#include <ios>
#include <iostream>
#include <vector>

#include <irtkImage.h>
#include <irtkFileToImage.h>
#include <irtkImageFunction.h>
#include <irtkResampling.h>

// Interpolation method
#define NEAREST_NEIGHBOR 0
#define LINEAR 1
#define BSPLINE 2
#define CSPLINE 3
#define SINC 4
#define SHAPE 5
#define GAUSSIAN 6

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

inline size_t index( size_t i, size_t j,
                     size_t shape0, size_t shape1 ) {
    return j + shape1*i;
}

inline size_t index( size_t i, size_t j, size_t k,
                     size_t shape0, size_t shape1, size_t shape2 ) {
    return k + shape2*( j + shape1*i );
}

inline size_t index( size_t i, size_t j, size_t k, size_t l,
                     size_t shape0, size_t shape1, size_t shape2, size_t shape3 ) {
    return l + shape3*( k + shape2*( j + shape1*i ) );
}

void Initialize();

int _get_header( char* filename,
                  double* pixelSize,
                  double* xAxis,
                  double* yAxis,
                  double* zAxis,
                  double* origin,
                  int* dim );

void put_attributes( irtkImageAttributes &attr,
                     double* pixelSize,
                     double* xAxis,
                     double* yAxis,
                     double* zAxis,
                     double* origin,
                     int* dim  );

void get_attributes( irtkImageAttributes &attr,
                     double* pixelSize,
                     double* xAxis,
                     double* yAxis,
                     double* zAxis,
                     double* origin,
                     int* dim );

template <class dtype>
void py2irtk( irtkGenericImage<dtype>& irtk_image,
              dtype* img,
              double* pixelSize,
              double* xAxis,
              double* yAxis,
              double* zAxis,
              double* origin,
              int* dim ) {
    
    irtkImageAttributes attr;

    put_attributes( attr,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );
        
    irtk_image.Initialize(attr);

    int n   = irtk_image.GetNumberOfVoxels();
    dtype* ptr = irtk_image.GetPointerToVoxels();
    for ( int i = 0; i < n; i++)
        ptr[i] = img[i];

    return;
}

template <class dtype>
void irtk2py( irtkGenericImage<dtype>& irtk_image,
              dtype* img,
              double* pixelSize,
              double* xAxis,
              double* yAxis,
              double* zAxis,
              double* origin,
              int* dim ) {
    
    irtkImageAttributes attr = irtk_image.GetImageAttributes();

    //int new_dim[4];
    get_attributes( attr,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );

    /* if ( dim[0] != new_dim[0] */
    /*      || dim[1] != new_dim[1] */
    /*      || dim[2] != new_dim[2] */
    /*      || dim[3] != new_dim[3] ) { */
    /*     std::cout << "the memory must be allocated in Python,\n" */
    /*               << "so dim must be known beforehand\n"; */
    /*     exit(1); */
    /* } */

    int n   = irtk_image.GetNumberOfVoxels();
    dtype* ptr = irtk_image.GetPointerToVoxels();
    for ( int i = 0; i < n; i++)
        img[i] = ptr[i];

    return;
}

template <class dtype>
void pyList2irtkVector( std::vector< irtkGenericImage<dtype> >& vec,
                        dtype* img,
                        double* pixelSize,
                        double* xAxis,
                        double* yAxis,
                        double* zAxis,
                        double* origin,
                        int* dim,
                        int n ) {

    // clean the vector first
    vec.clear();
    vec.reserve( n );

    size_t offset = 0;
    for ( size_t i = 0; i < n; i++ ) {
        irtkGenericImage<dtype> irtk_image;
        py2irtk<dtype>( irtk_image,
                        &img[offset],
                        &pixelSize[4*i],
                        &xAxis[3*i],
                        &yAxis[3*i],
                        &zAxis[3*i],
                        &origin[4*i],
                        &dim[4*i] );
        offset += dim[4*i]*dim[4*i+1]*dim[4*i+2]*dim[4*i+3];
        vec.push_back(irtk_image);
    }

}

template <class dtype>
void irtkVector2pyList( std::vector< irtkGenericImage<dtype> >& vec,
                        dtype* img,
                        double* pixelSize,
                        double* xAxis,
                        double* yAxis,
                        double* zAxis,
                        double* origin,
                        int* dim,
                        int& n ) {

    n = vec.size();

    size_t offset = 0;
    for ( size_t i = 0; i < n; i++ ) {
        py2irtk<dtype>( vec[i],
                        &img[offset],
                        &pixelSize[4*i],
                        &xAxis[3*i],
                        &yAxis[3*i],
                        &zAxis[3*i],
                        &origin[4*i],
                        &dim[4*i] );
        offset += dim[4*i]*dim[4*i+1]*dim[4*i+2]*dim[4*i+3];
    }

}

void _resample( float* img_in,
                double* pixelSize,
                double* xAxis,
                double* yAxis,
                double* zAxis,
                double* origin,
                int* dim,
                double* new_pixelSize,
                float* img_out,
                int interpolation_method=LINEAR,
                float gaussian_parameter=1.0 );

irtkMatrix py2matrix( int rows, int cols, double* data );

void _write_list( float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,
                 int n );

void _read_list( float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,
                 int& n );

void _transform_points( double* m, double* pts, size_t n );

void _points_to_image( unsigned char* img,
                       double* pixelSize,
                       double* xAxis,
                       double* yAxis,
                       double* zAxis,
                       double* origin,
                       int* dim,
                       double* pts,
                       size_t n );
    
#endif
