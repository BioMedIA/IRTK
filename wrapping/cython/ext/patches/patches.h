#include <cstdlib>
#include <iostream>
#include <cmath>

#define EPS 0.000001

typedef float pixel_t;

inline size_t index( size_t i, size_t j,
                     size_t shape0, size_t shape1 ) {
    return j + shape1*i;
}

inline size_t index( size_t i, size_t j, size_t k,
                     size_t shape0, size_t shape1, size_t shape2 ) {
    return k + shape2*( j + shape1*i );
}

inline size_t index( size_t i, size_t j, size_t k, size_t l,
                     size_t shape0, size_t shape1, size_t shape2, size_t shape3 ) {
    return l + shape3*( k + shape2*( j + shape1*i ) );
}

void _extract_patches2D( pixel_t* img,
                        int shape0, int shape1,
                        int* Y, int* X, size_t nb_points,
                        int rx, int ry,
                        pixel_t* patches );

void _project_gradients2D( double* grad,
                           int shape0, int shape1,
                           double* projection_matrix,
                           int nb_proj,
                           double* hist );

inline void integrate3D( double* sat,
                         int shape0, int shape1, int shape2,
                         int r0, int c0,
                         int r1, int c1,
                         double* res ) {
    for ( int k = 0; k < shape2; k++ )
        res[k] = sat[ index( r1, c1, k, shape0, shape1, shape2) ]
               + sat[ index( r0, c0, k, shape0, shape1, shape2) ]
               - sat[ index( r1, c0, k, shape0, shape1, shape2) ]
               - sat[ index( r0, c1, k, shape0, shape1, shape2) ];
}

void _extract_oriented_patches2D( double* img,
                                  int shape0, int shape1,
                                  double* hist,
                                  double* projection_matrix,
                                  int nb_proj,
                                  int* X, int* Y, size_t nb_points,
                                  int r,
                                  double* patches );
