#include "irtk2cython.h"

void Initialize() {
    // turn off the synchronization of iostream objects and cstdio streams
    // for increased speed
    std::ios_base::sync_with_stdio(false);
}

int _get_header( char* filename,
                  double* pixelSize,
                  double* xAxis,
                  double* yAxis,
                  double* zAxis,
                  double* origin,
                  int* dim ) {
    try {
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
    catch ( irtkException &e ) {
        return -1;
    }
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
    
    irtkResampling<float> resampling( new_pixelSize[0],
                                      new_pixelSize[1],
                                      new_pixelSize[2] );
    resampling.SetInput ( &irtk_image );
    resampling.SetOutput( &irtk_image );
    resampling.SetInterpolator( interpolator );
    resampling.Run();

    irtk2py<float>( irtk_image,
                    img_out,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );
}

void _gaussianBlurring( float* img_in,
                        double* pixelSize,
                        double* xAxis,
                        double* yAxis,
                        double* zAxis,
                        double* origin,
                        int* dim,
                        float* img_out,
                        double sigma ) {

    irtkGenericImage<float> irtk_image;
    py2irtk<float>( irtk_image,
                    img_in,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );

    irtkGaussianBlurring<float> gb(sigma);
    gb.SetInput(&irtk_image);
    gb.SetOutput(&irtk_image);
    gb.Run();

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

// test function
void _write_list( float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,
                 int n ) {
    std::vector< irtkGenericImage<float> > vec;
    pyList2irtkVector( vec,
                       img,
                       pixelSize,
                       xAxis,
                       yAxis,
                       zAxis,
                       origin,
                       dim,
                       n );

    
    for ( int i = 0; i < n; i++ ) {
        std::stringstream ss;
        ss << "file" << i << ".nii";
        vec[i].Write( ss.str().c_str() );
    }
    
}

// test function
void _read_list( float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,
                 int& n ) {
    // we need to know the number of pixels beforehand,
    // or a maximal number of pixels for Python to allocate memory...
    std::vector< irtkGenericImage<float> > vec;
    for ( int i = 0; i < 3; i++ ) {
        irtkGenericImage<float> irtk_image;
        std::stringstream ss;
        ss << "file" << i << ".nii";
        irtk_image.Read( ss.str().c_str() );
        vec.push_back( irtk_image );
    }
    irtkVector2pyList( vec,
                       img,
                       pixelSize,
                       xAxis,
                       yAxis,
                       zAxis,
                       origin,
                       dim,
                       n );

}

void _transform_points( double* m, double* pts, size_t n ) {
    double tmp_pt[3];
    for ( size_t i = 0; i < n; i++ ) {
        tmp_pt[0] = m[index(0,0,4,4)] * pts[index(i,0,n,3)] + m[index(0,1,4,4)] * pts[index(i,1,n,3)] + m[index(0,2,4,4)] * pts[index(i,2,n,3)] + m[index(0,3,4,4)];
        tmp_pt[1] = m[index(1,0,4,4)] * pts[index(i,0,n,3)] + m[index(1,1,4,4)] * pts[index(i,1,n,3)] + m[index(1,2,4,4)] * pts[index(i,2,n,3)] + m[index(1,3,4,4)];
        tmp_pt[2] = m[index(2,0,4,4)] * pts[index(i,0,n,3)] + m[index(2,1,4,4)] * pts[index(i,1,n,3)] + m[index(2,2,4,4)] * pts[index(i,2,n,3)] + m[index(2,3,4,4)];
        pts[index(i,0,n,3)] = tmp_pt[0];
        pts[index(i,1,n,3)] = tmp_pt[1];
        pts[index(i,2,n,3)] = tmp_pt[2];
    }        
}

void _points_to_image( unsigned char* img,
                       double* pixelSize,
                       double* xAxis,
                       double* yAxis,
                       double* zAxis,
                       double* origin,
                       int* dim,
                       double* pts,
                       size_t n ) {

    irtkGenericImage<unsigned char> irtk_image;
    py2irtk<unsigned char>( irtk_image,
                    img,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );

    double pt[3];
    for ( size_t i = 0; i < n; i++ ) {
        pt[0] = pts[index(i,0,n,3)];
        pt[1] = pts[index(i,1,n,3)];
        pt[2] = pts[index(i,2,n,3)];
        irtk_image.WorldToImage( pt[0], pt[1], pt[2] );
        if ( pt[0] >= 0 && pt[0] < irtk_image.GetX()
             && pt[1] >= 0 && pt[1] < irtk_image.GetY()
             && pt[2] >= 0 && pt[2] < irtk_image.GetZ() )
            irtk_image( pt[0], pt[1], pt[2] ) = 1;
    }

    irtk2py<unsigned char>( irtk_image,
                    img,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );    
}

void _read_points( char *filename,
                   std::vector< std::vector<double> > &points ) {
    irtkPointSet irtk_points;
    irtk_points.ReadVTK( filename );
    for ( int i = 0; i < irtk_points.Size(); i++ ) {
        irtkPoint p;
        p = irtk_points(i);
        std::vector<double> pt(3);
        pt[0] = p._x; pt[1] = p._y; pt[2] = p._z;
        points.push_back( pt );
    }
}

void _write_points( char *filename,
                    std::vector< std::vector<double> > &points ) {
    irtkPointSet irtk_points;
    for ( int i = 0; i < irtk_points.Size(); i++ ) {
        irtk_points.Add( irtkPoint( points[i][0],
                                    points[i][1],
                                    points[i][2] ) );
    }
    irtk_points.WriteVTK( filename );
}

void _orientation( double* pixelSize,
                   double* xAxis,
                   double* yAxis,
                   double* zAxis,
                   double* origin,
                   int* dim,
                   int &i, int &j, int &k) {

    irtkGenericImage<uchar> irtk_image;
    irtkImageAttributes attr;

    put_attributes( attr,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );
        
    irtk_image.Initialize(attr);

    irtk_image.Orientation( i, j, k );
}
