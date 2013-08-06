#include "patches.h"

void _extract_patches2D( pixel_t* img,
                        int shape0, int shape1,
                        int* X, int* Y, size_t nb_points,
                        int rx, int ry,
                        pixel_t* patches ) {
                    std::cout 
                              <<  nb_points<< " "
                              <<  2*ry+1<< " "
                              <<2*rx+1<< "\n";      
    for ( size_t i = 0; i < nb_points; i++ )
        for ( int x = X[i] - rx; x <= X[i] + rx; x++ )
            for ( int y = Y[i] - ry; y <= Y[i] + ry; y++ ) {
               
                if ( 0 <= x && x < shape1
                     && 0 <= y && y < shape0 ) {
                    // std::cout << i << " "
                    //           << y-(Y[i]-ry)<< " "
                    //           <<   x-(X[i]-rx)<< " "
                    //           <<  nb_points<< " "
                    //           <<  2*ry+1<< " "
                    //           <<2*rx+1<< "\n";
                    
                    patches[index(i,
                                  y-(Y[i]-ry),
                                  x-(X[i]-rx),
                                  nb_points,
                                  2*ry+1,
                                  2*rx+1)] = img[index(y,x,shape0,shape1)];
                }
            }
}

void _reconstruct_from_patches2D() {

}

void _project_gradients2D( double* grad,
                           int shape0, int shape1,
                           double* projection_matrix,
                           int nb_proj,
                           double* hist ) {    
    for ( int i = 0; i < shape0; i++ )
        for ( int j = 0; j < shape1; j++ )
            for ( int a = 0; a < nb_proj; a++ )
                hist[index(i,j,a, shape0, shape1, nb_proj)]
                    = projection_matrix[index(a,0,nb_proj,2)]
                    * grad[index(i,j,0,shape0,shape1,2)]
                    + projection_matrix[index(a,1,nb_proj,2)]
                    * grad[index(i,j,1,shape0,shape1,2)];
}

int argmax( double* a, int shape0 ) {
    int i_max = 0;
    double m = a[0];
    for ( int i = 1; i < shape0; i++ )
        if (a[i] > m) {
            m = a[i];
            i_max = i;
        }
    return i_max;
}

void orthogonal_vector( double* a,
                        double* res ) {
    res[0] = -a[1];
    res[1] = a[0];
}

void get_orientation2D( double* hist,
                        int shape0, int shape1,
                        int nb_proj,
                        double* projection_matrix,
                        int x, int y, int r,
                        double* u, double* v) {

    int x0 = x - r;
    int y0 = y - r;

    int x0_end = x + r;
    int y0_end = y + r;
    
    double* local_hist = new double[nb_proj];
    integrate3D( hist,
                 shape0, shape1, nb_proj,
                 y0, x0,
                 y0_end, x0_end,
                 local_hist );
    
    int direction = argmax( local_hist, nb_proj );

    for ( int d = 0; d < 2; d++ )
        u[d] = projection_matrix[index(direction,d,nb_proj,2)];
       
    orthogonal_vector( u, v );

    delete local_hist;

}

void extract_patch2D( double* img,
                      int shape0, int shape1,
                      double* u, double* v,
                      int x, int y, int r,
                      double* patch ) {

    int p_shape0 = 2*r + 1;
    int p_shape1 = 2*r + 1;

    double offset_x = -r - 1;
    double offset_y = -r - 1;

    double start_y = y
        + offset_y*u[0]
        + offset_x*v[0];     
    double start_x = x
        + offset_y*u[1]
        + offset_x*v[1];

    for ( int i = 0; i < p_shape0; i++ ) {
        start_y += u[0];
        start_x += u[1];

        double pixel_y = start_y;
        double pixel_x = start_x;
        for ( int j = 0; j < p_shape1; j++ ) {
            pixel_y += v[0];
            pixel_x += v[1];

            // i0 = std::min(std::max(round(pixel_y), 0), (int)shape0-1);
            // j0 = std::min(std::max(round(pixel_x), 0), (int)shape1-1);
                
            int i0 = pixel_y;
            int j0 = pixel_x;
                
            if ( i0 < 0 || i0 >= shape0
                 || j0 < 0 || j0 >= shape1 ) {
                patch[index(i,j,p_shape0,p_shape1)] = 0;
            }
            else {
                // Nearest neighbour interpolation
                /*
                patch[index(i,j,p_shape0,p_shape1)]
                    = img[index(i0,j0,shape0,shape1)];
                */

                // Linear interpolation
                double u2 = pixel_y - i0;
                double v2 = pixel_x - j0;
                double u1 = 1 - u2;
                double v1 = 1 - v2;

                double val00, val01, val10, val11;
                val00 = img[index(i0,j0,shape0,shape1)];
                if ( i0 == shape0 - 1) {
                    val10 = 0;
                    u1 = 1; u2 = 0;
                }
                else {
                    val10 = img[index(i0+1,j0,shape0,shape1)];
                }
                if ( j0 == shape1 - 1) {
                    val01 = 0;
                    v1 = 1; v2 = 0;
                }
                else {
                    val01 = img[index(i0,j0+1,shape0,shape1)];
                }
                if ( ( i0 == shape0 - 1) || ( j0 == shape1 - 1) ) {
                    val11 = 0;
                }
                else {
                    val11 = img[index(i0+1,j0+1,shape0,shape1)];
                }
                    
                patch[index(i,j,p_shape0,p_shape1)] =
                    u1*v1*val00 +
                    u2*v1*val10 +
                    u1*v2*val01 +
                    u2*v2*val11;
            }
        }
    }
}

void _extract_oriented_patches2D( double* img,
                                  int shape0, int shape1,
                                  double* hist,
                                  double* projection_matrix,
                                  int nb_proj,
                                  int* X, int* Y, size_t nb_points,
                                  int r,
                                  double* patches ) {

    int patch_voxels = (2*r+1)*(2*r+1);
    double u[2];
    double v[2];
    double* patch = new double[patch_voxels];
    for ( int i = 0; i < nb_points; i++ ) {

        if ( Y[i] < r || X[i] < r
             || Y[i] >= shape0-r
             || X[i] >= shape1-r )
            continue;

        get_orientation2D( hist,
                           shape0+1, shape1+1,
                           nb_proj,
                           projection_matrix,
                           X[i], Y[i], r,
                           u, v);

        extract_patch2D( img,
                         shape0, shape1,
                         u, v,
                         X[i], Y[i], r,
                         patch );

        // store patch
        for ( int d = 0; d < patch_voxels; d++ )
            patches[index(i,d,nb_points,patch_voxels)] = patch[d];

    }
    
    delete patch;
}

/*
double descriptor2D( double* patch,
                     double* descr,
                     int r,
                     int NB_SUBREGIONS,
                     int SUBREGION_SIZE ) {

    int descr_length = 6*NB_SUBREGIONS*NB_SUBREGIONS*NB_SUBREGIONS;
    int PATCH_SIZE = 2*r + 1;
    
    for ( int i = 0; i < NB_SUBREGIONS; i++ )
        for ( int j = 0; j < NB_SUBREGIONS; j++ ) {
            double dy = 0;
            double dx = 0;
            for ( int y = i*SUBREGION_SIZE; y < (i+1)*SUBREGION_SIZE; y++ )
                for ( int x = j*SUBREGION_SIZE; x < (j+1)*SUBREGION_SIZE; x++ ) {

                    dy = patch[index(y,x,PATCH_SIZE,PATCH_SIZE)]
                        - patch[index(y+1,x,PATCH_SIZE,PATCH_SIZE)]
                        + patch[index(y,x+1,PATCH_SIZE,PATCH_SIZE)]
                        - patch[index(y+1,x+1,PATCH_SIZE,PATCH_SIZE)];

                    descr[index(i,j,0,NB_SUBREGIONS,NB_SUBREGIONS,4)] += dy;
                    descr[index(i,j,1,NB_SUBREGIONS,NB_SUBREGIONS,4)] += fabs(dy);
                    
                    dx = patch[index(y,x,PATCH_SIZE,PATCH_SIZE)]
                        - patch[index(y,x+1,PATCH_SIZE,PATCH_SIZE)]
                        + patch[index(y+1,x,PATCH_SIZE,PATCH_SIZE)]
                        - patch[index(y+1,x+1,PATCH_SIZE,PATCH_SIZE)];

                    descr[index(i,j,2,NB_SUBREGIONS,NB_SUBREGIONS,4)] += dx;
                    descr[index(i,j,3,NB_SUBREGIONS,NB_SUBREGIONS,4)] += fabs(dx);
                }   
        }

    float norm = 0;
    for ( int d = 0; d < descr_length; d++ )
        norm += descr[d]*descr[d];
    if (norm > 0) {
        norm = 255 * 1.0/sqrt(norm + EPS);
        for ( int d = 0; d < descr_length; d++ )
            descr[d] *= norm;
    }

    return norm;
}
*/

