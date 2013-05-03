#ifndef IRTK2CYTHON_H
#define IRTK2CYTHON_H

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

    get_attributes( attr,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );

    int n   = irtk_image.GetNumberOfVoxels();
    dtype* ptr = irtk_image.GetPointerToVoxels();
    for ( int i = 0; i < n; i++)
        img[i] = ptr[i];

    return;
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


#endif
