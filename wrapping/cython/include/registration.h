#ifndef REGISTRATION_H
#define REGISTRATION_H


#include <irtkTransformation.h>
#include <irtkRegistration.h>

#include "irtk2cython.h"
void rigid2py( irtkRigidTransformation &transform,
               double &tx,
               double &ty,
               double &tz,
               double &rx,
               double &ry,
               double &rz,
               bool invert=false );

void py2rigid( irtkRigidTransformation &transform,
               double tx,
               double ty,
               double tz,
               double rx,
               double ry,
               double rz,
               bool invert=false );

void pyList2rigidVector( std::vector< irtkRigidTransformation > &vec,
                         double* tx,
                         double* ty,
                         double* tz,
                         double* rx,
                         double* ry,
                         double* rz,
                         int n,
                         bool invert=false );

void rigidVector2pyList( std::vector< irtkRigidTransformation > &vec,
                         double* tx,
                         double* ty,
                         double* tz,
                         double* rx,
                         double* ry,
                         double* rz,
                         bool invert=false );

void _read_rigid( char* filename,
                  double &tx,
                  double &ty,
                  double &tz,
                  double &rx,
                  double &ry,
                  double &rz );

void _write_rigid( char* filename,
                   double tx,
                   double ty,
                   double tz,
                   double rx,
                   double ry,
                   double rz );

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
                       int interpolation_method=LINEAR,
                       float gaussian_parameter=1.0 );

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
                          double &rz );

double _registration_rigid_points( double* source_points,
                                   double* target_points,
                                   int n,
                                   double &tx,
                                   double &ty,
                                   double &tz,
                                   double &rx,
                                   double &ry,
                                   double &rz );

#endif
