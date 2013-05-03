#include "irtk2cython.h"

int _get_header( char* filename,
                  double* pixelSize,
                  double* xAxis,
                  double* yAxis,
                  double* zAxis,
                  double* origin,
                  int* dim ) {
    // irtkGreyImage image(filename);
    // image.GetPixelSize(&pixelSize[0], &pixelSize[1], &pixelSize[2]);
    // image.GetOrientation(xAxis, yAxis, zAxis);
    // image.GetOrigin(origin[0], origin[1], origin[2]);
    // dim[0] = image.GetX();
    // dim[1] = image.GetY();
    // dim[2] = image.GetZ();
    // dim[3] = image.GetT();

    irtkBaseImage *image;
    irtkFileToImage *reader = irtkFileToImage::New(filename);
    image = reader->GetOutput();
    image->GetPixelSize(&pixelSize[0], &pixelSize[1], &pixelSize[2], &pixelSize[3]);
    image->GetOrientation(xAxis, yAxis, zAxis);
    image->GetOrigin(origin[0], origin[1], origin[2], origin[3]);
    dim[0] = image->GetX();
    dim[1] = image->GetY();
    dim[2] = image->GetZ();
    dim[3] = image->GetT();

    int dtype = reader->GetDataType();

    delete reader;
    delete image;

    return dtype;
}

void put_attributes( irtkImageAttributes &attr,
                     double* pixelSize,
                     double* xAxis,
                     double* yAxis,
                     double* zAxis,
                     double* origin,
                     int* dim  ) {
    // dim
    attr._x = dim[0];
    attr._y = dim[1];
    attr._z = dim[2];
    attr._t = dim[3];

    // voxel size
    attr._dx = pixelSize[0];
    attr._dy = pixelSize[1];
    attr._dz = pixelSize[2];
    attr._dt = pixelSize[3];
    
    // origin
    attr._xorigin = origin[0];
    attr._yorigin = origin[1];
    attr._zorigin = origin[2];
    attr._torigin = origin[3];

    // x-axis
    attr._xaxis[0] = xAxis[0];
    attr._xaxis[1] = xAxis[1];
    attr._xaxis[2] = xAxis[2];

    // y-axis
    attr._yaxis[0] = yAxis[0];
    attr._yaxis[1] = yAxis[1];
    attr._yaxis[2] = yAxis[2];

    // z-axis
    attr._zaxis[0] = zAxis[0];
    attr._zaxis[1] = zAxis[1];
    attr._zaxis[2] = zAxis[2];
}

void get_attributes( irtkImageAttributes &attr,
                     double* pixelSize,
                     double* xAxis,
                     double* yAxis,
                     double* zAxis,
                     double* origin,
                     int* dim ) {
    // dim
    dim[0] = attr._x;
    dim[1] = attr._y;
    dim[2] = attr._z;
    dim[3] = attr._t;

    // voxel size
    pixelSize[0] = attr._dx;
    pixelSize[1] = attr._dy;
    pixelSize[2] = attr._dz;
    pixelSize[3] = attr._dt;
    
    // origin
    origin[0] = attr._xorigin;
    origin[1] = attr._yorigin;
    origin[2] = attr._zorigin;
    origin[3] = attr._torigin;

    // x-axis
    xAxis[0] = attr._xaxis[0];
    xAxis[1] = attr._xaxis[1];
    xAxis[2] = attr._xaxis[2];

    // y-axis
    yAxis[0] = attr._yaxis[0];
    yAxis[1] = attr._yaxis[1];
    yAxis[2] = attr._yaxis[2];

    // z-axis
    zAxis[0] = attr._zaxis[0];
    zAxis[1] = attr._zaxis[1];
    zAxis[2] = attr._zaxis[2];
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
                int interpolation_method,
                float gaussian_parameter ) {
    
    irtkGenericImage<float> irtk_image;// = new irtkGenericImage<float>;
    irtkImageFunction *interpolator = NULL;

    switch (interpolation_method) {
          
    case NEAREST_NEIGHBOR:
        { interpolator = new irtkNearestNeighborInterpolateImageFunction; }
        break;

    case LINEAR:
      { interpolator = new irtkLinearInterpolateImageFunction; }
      break;

    case BSPLINE:
      { interpolator = new irtkBSplineInterpolateImageFunction; }
      break;

    case CSPLINE:
        { interpolator = new irtkCSplineInterpolateImageFunction; }
      break;

    case SINC:
        { interpolator = new irtkSincInterpolateImageFunction; }
      break;

    case SHAPE:
        { interpolator = new irtkShapeBasedInterpolateImageFunction; }
        break;

    case GAUSSIAN:
        { interpolator = new irtkGaussianInterpolateImageFunction(gaussian_parameter); }
        break;

      default:
      cout << "Unknown interpolation method" << endl;
  }
    
    py2irtk<float>( irtk_image,
                    img_in,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );

    irtk_image.Print();
    
    irtkResampling<float> resampling( new_pixelSize[0],
                                      new_pixelSize[1],
                                      new_pixelSize[2] );
    resampling.SetInput ( &irtk_image );
    resampling.SetOutput( &irtk_image );
    resampling.SetInterpolator( interpolator );
    resampling.Run();

    irtk_image.Print();

    irtk2py<float>( irtk_image,
                    img_out,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );
}

irtkMatrix py2matrix( int rows, int cols, double* data ) {
    irtkMatrix matrix(rows, cols);
    int i, j;
    for ( i = 0; i < rows; i++ )
        for ( j = 0; j < cols; j++ )
            matrix(i,j) = data[index(i,j,rows,cols)];

    return matrix;
}

