#include <stdio.h>
#include <iostream>
#include <cmath>
#include "GCoptimization.h"

typedef float pixel_t;

typedef int LabelID;
typedef int SiteID;
typedef float EnergyTermType;
typedef float EnergyType;

inline size_t index( size_t i, size_t j,
                     size_t shape0, size_t shape1 ) {
    return j + shape1*i;
}

inline size_t index( size_t i, size_t j, size_t k,
                     size_t shape0, size_t shape1, size_t shape2 ) {
    return k + shape2*( j + shape1*i );
}

void _crf( pixel_t* img,
           int shape0, int shape1,
           double std,
           LabelID* labels,
           double* proba,
           double lambda,
           int degree );

void _crf3D( pixel_t* img,
          int shape0, int shape1, int shape2,
          double std,
          LabelID* labels,
          double* proba,
             double lambda );

