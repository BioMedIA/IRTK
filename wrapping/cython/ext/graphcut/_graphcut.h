#include <stdio.h>
#include <iostream>
#include <cmath>
#include "graph.h"

typedef double voxel_t;

typedef double captype;
typedef double tcaptype;
typedef double flowtype;

inline size_t index( size_t i, size_t j, size_t k,
                     size_t shape0, size_t shape1, size_t shape2 ) {
    return k + shape2*( j + shape1*i );
}

void _graphcut( voxel_t* img,
               int shape0, int shape1, int shape2,
               double img_std,
               unsigned char* mask,
               double* raw_spacing,
               unsigned char* seg );
