/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCommon.h>

#include <nr.h>

double median(float q, float p0, float p1){
    float m;
    m = p0>p1?p0:p1;
    m = q > m? m:q;
    return m;
}

/// weight must be normalized to sum = 1 before entering here.
/// implemented fast weighted median selector algorithm from Y. Li and S. Osher.
double weightedmedian(int index, double lambda, double f, float *neighbors, float *weight)
{
    int i,pos;
    double p0,p1,csum;

    // sort the value according to neighbors
    sort2(index-1,neighbors,weight);
    
    csum = 1;

    if(lambda > 0){

        p0 = f+(2.0*csum/lambda);
        // find new median
        for(i = 1; i < index; i++){
             csum -= 2.0*weight[i];
             p1 = f+(2.0*csum/lambda);
             p0 = median(p0,p1,neighbors[i]);
        }

        return p0;
    }else{
        // 1/lambda = infinite, use a simple version.
        pos = 0;

        for(i = 1; i <= index; i++){
            if(csum > 0) pos++;
            if(csum < 0) pos--;
            if(i != index){
                csum = csum - 2.0*weight[i];
            }
        }
        return neighbors[(index+pos)/2];
    }
}