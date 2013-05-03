#ifndef VOXELLISE_H
#define VOXELLISE_H

#include "irtk2cython.h"

#include <irtkImage.h>

#include <irtkGaussianBlurring.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkStripper.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencil.h>
#include <vtkImageData.h>
#include <vtkTriangleFilter.h>
#include <vtkStructuredPoints.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>

void voxellise( vtkPolyData *poly, // input mesh (must be closed)
                irtkGreyImage &image,
                double value=1 );

void create_polydata( double* points,
                      int npoints,
                      int* triangles,
                      int ntriangles,
                      vtkPolyData *poly );

void _voxellise( double* points,
                 int npoints,
                 int* triangles,
                 int ntriangles,
                 short* img, // irtkGreyImage
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim );
#endif
