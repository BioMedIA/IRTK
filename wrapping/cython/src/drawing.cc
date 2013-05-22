#include "drawing.h"

void _drawSphere( unsigned char* img,
                  int shape0,
                  int shape1,
                  int shape2,
                  int x0,
                  int y0,
                  int z0,
                  int rx,
                  int ry,
                  int rz ) {
    int rx2 = rx*rx;
    int ry2 = ry*ry;
    int rz2 = rz*rz;
    int x, y, z;
    for ( int dz = -rz; dz <= rz; dz++ )
        for ( int dy = -ry; dy <= ry; dy++ )
            for ( int dx = -rx; dx <= rx; dx++ ) {                
                x = x0 + dx;
                y = y0 + dy;
                z = z0 + dz;
                if ( z >= 0 && z < shape0
                     && y >= 0 && y < shape1
                     && x >= 0 && x < shape2 )
                    if ( float(dx*dx)/rx2
                         + float(dy*dy)/ry2
                         + float(dz*dz)/rz2 <= 1 )
                        img[index(z,y,x,shape0,shape1,shape2)] = 1;                 
            }                    
                    
                    // if ( x*x + y*y + z*z <= rx2 ) {
                    //     cout << "setting to 1\n";
                    //     img[index(z,y,x,shape0,shape1,shape2)] = 1;
                    // }
}
