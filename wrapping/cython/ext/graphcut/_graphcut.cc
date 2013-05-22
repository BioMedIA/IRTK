#include "_graphcut.h"

inline double _pow( double x, int a ) { return x*x; }

void _graphcut( voxel_t* img,
               int shape0, int shape1, int shape2,
               double img_std,
               unsigned char* mask,
               double* raw_spacing,
               unsigned char* seg ) {

    typedef Graph<captype,tcaptype,flowtype> GraphType;

    size_t nb_voxels = shape0*shape1*shape2;

    std::cout << "will create graph... " << nb_voxels  << "\n";
    GraphType *G = new GraphType( /*estimated # of nodes*/ nb_voxels,
                                  /*estimated #	of edges*/ nb_voxels*(26+1) );

    G->add_node( nb_voxels );

    std::cout << "Image std:" << img_std << "\n";
    double img_std2 = 2 * _pow( img_std, 2);
    
    std::cout << "building graph...\n";

    int i,j,k,a,b,c;
    size_t id, id0;
    int i0, j0, k0;
    double dist, w;
    for ( i = 0; i < shape0; i++ )
        for ( j = 0; j < shape1; j++ )
            for ( k = 0; k < shape2; k++ )
                for ( a = -1; a <= 1; a++ )
                    for ( b = -1; b <= 1; b++ )
                        for ( c = -1; c <= 1; c++ ) {
                            if ( abs(a)+abs(b)+abs(c) == 0 )
                                continue;
                            i0 = i+a;
                            j0 = j+b;
                            k0 = k+c;
                            if ( 0 <= i0  && i0 < shape0
                                 && 0 <= j0 && j0 < shape1
                                 && 0 <= k0 && k0 < shape2 ) {
                                id = index(i,j,k,shape0,shape1,shape2);
                                id0 = index(i0,j0,k0,shape0,shape1,shape2);
                                dist = sqrt( _pow( (double)(a)*raw_spacing[2], 2 )
                                             + _pow( (double)(b)*raw_spacing[1], 2 )
                                             + _pow( (double)(c)*raw_spacing[0], 2 )
                                             );
                                if ( img[id] < img[id0] )
                                    w = 1.0/dist;
                                else
                                    w = exp( - _pow(img[id] - img[id0], 2)
                                             / img_std2 )
                                        / dist;
                                //std::cout << w << " " << dist <<"\n";
                                G->add_edge( id,
                                             id0,
                                             w,
                                             0);
                                    }
                        }

    std::cout << "linking to source and sink...\n";
    
    for ( i = 0; i < shape0; i++ )
        for ( j = 0; j < shape1; j++ )
            for ( k = 0; k < shape2; k++ ) {
                id = index(i,j,k,shape0,shape1,shape2);
                if ( mask[id] == 1 )
                    G->add_tweights(id,1000,0);
                else if ( mask[id] == 2 ) 
                    G->add_tweights(id,0,1000);
            }

    std::cout << "computing maxflow...\n";
    std::cout << G->maxflow() << "\n";

    std::cout <<  "transcription of segmentation...\n";
    
    for ( i = 0; i < shape0; i++ )
        for ( j = 0; j < shape1; j++ )
            for ( k = 0; k < shape2; k++ ) {
                id = index(i,j,k,shape0,shape1,shape2);
                if ( G->what_segment(id) == 0 )
                    seg[id] = 1;
            }

    return;
}
