#include "crf.h" 

void _crf( pixel_t* img,
           double* pixelSize,
           double* xAxis,
           double* yAxis,
           double* zAxis,
           double* origin,
           int* dim,
           LabelID* labels,
           double* proba,
           double l ) {

    irtkGenericImage<pixel_t> irtk_image;
    py2irtk<pixel_t>( irtk_image,
                      img,
                      pixelSize,
                      xAxis,
                      yAxis,
                      zAxis,
                      origin,
                      dim );

    irtkGenericImage<LabelID> irtk_labels;
    py2irtk<LabelID>( irtk_labels,
                      labels,
                      pixelSize,
                      xAxis,
                      yAxis,
                      zAxis,
                      origin,
                      dim );

    irtkGenericImage<double> irtk_proba;
    py2irtk<double>( irtk_proba,
                     proba,
                     pixelSize,
                     xAxis,
                     yAxis,
                     zAxis,
                     origin,
                     dim );

    irtkCRF crf( irtk_image,
                 irtk_labels,
                 irtk_proba );
    crf.SetLambda( l );
    crf.Run();

    irtk2py<LabelID>( irtk_labels,
                      labels,
                      pixelSize,
                      xAxis,
                      yAxis,
                      zAxis,
                      origin,
                      dim );
}
