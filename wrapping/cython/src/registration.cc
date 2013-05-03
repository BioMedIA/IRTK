#include "registration.h"

void rigid2py( irtkRigidTransformation &transform,
               double &tx,
               double &ty,
               double &tz,
               double &rx,
               double &ry,
               double &rz ) {
    tx = transform.GetTranslationX();
    ty = transform.GetTranslationY();
    tz = transform.GetTranslationZ();
    rx = transform.GetRotationX();
    ry = transform.GetRotationY();
    rz = transform.GetRotationZ();
}

void py2rigid( irtkRigidTransformation &transform,
               double tx,
               double ty,
               double tz,
               double rx,
               double ry,
               double rz ) {
    transform.PutTranslationX( tx );
    transform.PutTranslationY( ty );
    transform.PutTranslationZ( tz );
    transform.PutRotationX( rx );
    transform.PutRotationY( ry );
    transform.PutRotationZ( rz );
}

void _read_rigid( char* filename,
                  double &tx,
                  double &ty,
                  double &tz,
                  double &rx,
                  double &ry,
                  double &rz ) {
    irtkRigidTransformation transform;
    transform.irtkTransformation::Read( filename );
    rigid2py( transform,
              tx,ty,tz,
              rx ,ry, rz );
}

void _write_rigid( char* filename,
                   double tx,
                   double ty,
                   double tz,
                   double rx,
                   double ry,
                   double rz ) {
    irtkRigidTransformation transform;
    py2rigid( transform,
              tx,ty,tz,
              rx ,ry, rz );
    transform.irtkTransformation::Write( filename );
}

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
                       float gaussian_parameter ) {

    // transformation
    irtkRigidTransformation transform;
    py2rigid( transform,
              tx,ty,tz,
              rx ,ry, rz );

    // source
    irtkGenericImage<float> source;
    py2irtk<float>( source,
                    source_img,
                    source_pixelSize,
                    source_xAxis,
                    source_yAxis,
                    source_zAxis,
                    source_origin,
                    source_dim );

    // target
    irtkGenericImage<float> target;   
    irtkImageAttributes attr;
    put_attributes( attr,
                    target_pixelSize,
                    target_xAxis,
                    target_yAxis,
                    target_zAxis,
                    target_origin,
                    target_dim );
    target.Initialize( attr );

    // interpolator
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

    // Create image transformation
    irtkImageTransformation imagetransformation;

    imagetransformation.SetInput( &source, &transform );
    imagetransformation.SetOutput( &target );
    imagetransformation.PutInterpolator( interpolator );

    // padding
    int target_padding, source_padding;
    source_padding = 0;
    target_padding = MIN_GREY;
    imagetransformation.PutTargetPaddingValue( target_padding );
    imagetransformation.PutSourcePaddingValue( source_padding );

    // inverse transformation
    // imagetransformation.InvertOn();
    
    // Transform image
    imagetransformation.Run();

    irtk2py<float>( target,
                    target_img,
                    target_pixelSize,
                    target_xAxis,
                    target_yAxis,
                    target_zAxis,
                    target_origin,
                    target_dim );

    delete interpolator;
}

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
                          double &rz ) {
    
    // source
    /** Second input image. This image is denoted as source image. The goal of
     *  the registration is to find the transformation which maps the source
     *  image into the coordinate system of the target image.
     */
    irtkGenericImage<short> source;
    py2irtk<short>( source,
                    source_img,
                    source_pixelSize,
                    source_xAxis,
                    source_yAxis,
                    source_zAxis,
                    source_origin,
                    source_dim );
    

    // target
    /** First set of input image. This image is denoted as target image and its
     *  coordinate system defines the frame of reference for the registration.
     */
    irtkGenericImage<short> target;
    py2irtk<short>( target,
                    target_img,
                    target_pixelSize,
                    target_xAxis,
                    target_yAxis,
                    target_zAxis,
                    target_origin,
                    target_dim );     
    
    // Create transformation
    irtkRigidTransformation transformation;

    // Initialize transformation
    py2rigid( transformation,
              tx,ty,tz,
              rx ,ry, rz );

    transformation.Invert(); // inverts only the matrix
    transformation.UpdateParameter();
    
    // Create registration
    // The goal of the registration is to find the transformation which maps the
    // source image into the coordinate system of the target image whereas
    // irtkImageTransformation expect the inverse transformations, hence the
    // calls to Invert().
    irtkImageRigidRegistrationWithPadding registration;
    
    registration.SetInput( &source, &target );
    registration.SetOutput( &transformation );

    // Make an initial Guess for the parameters.
    //registration.GuessParameter();
    registration.GuessParameterSliceToVolume();
    
    // TODO: Overrride with any the user has set.
    // use parameter file?
    
    // Run registration filter
    registration.Run();

    transformation.Invert(); // inverts only the matrix
    transformation.UpdateParameter();

    // We return the transformation mapping locations in the target image to
    // locations in the source image (this is the input expected by
    // irtkImageTransformation). 

    rigid2py( transformation,
              tx,ty,tz,
              rx ,ry, rz );
    
}

