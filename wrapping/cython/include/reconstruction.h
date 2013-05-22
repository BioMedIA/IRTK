#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "irtk2cython.h"
#include "registration.h"

#include <irtkReconstruction.h>

void _reconstruct(
                 // input stacks or slices
                 float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,

                 // number of stacks
                 int n,

                 // stack ids: which stack each slice
                 // comes from
                 int* _stack_ids,

                 // number of reconstruction iterations to run
                 int iterations,

                 // initial transformations
                 double* tx,
                 double* ty,
                 double* tz,
                 double* rx,
                 double* ry,
                 double* rz,

                 // slice thickness
                 double* _thickness,

                 // mask (header same as template)
                 float* mask_img,
                 
                 // output: reconstructed image
                 float* reconstructed_img,
                 double* reconstructed_pixelSize,
                 double* reconstructed_xAxis,
                 double* reconstructed_yAxis,
                 double* reconstructed_zAxis,
                 double* reconstructed_origin,
                 int* reconstructed_dim );

#endif
