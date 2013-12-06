#include "irtk2cython.h"
#include "irtkCRF.h"

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
           double sigmaZ );
