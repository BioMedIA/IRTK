/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkCommon.h>

#include <algorithm>

double median(float q, float p0, float p1){
    float m;
    m = p0>p1?p0:p1;
    m = q > m? m:q;
    return m;
}

void sort2(int index, float *array1, float *array2) {
    int* indices = new int[index];
    float* array1_store = new float[index];
    float* array2_store = new float[index];

    for (int i = 0; i < index; i++) {
       // create array with indices to be sorted
       indices[i] = i;
       // keep a back up of the original arrays
       array1_store[i] = array1[i];
       array2_store[i] = array2[i];
    }

    // sort the first array and keep the index order
    std::sort(indices, indices + index, sort_indices(array1));

    // update the arrays
    for (int i = 0; i < index; i++) {
       array1[i+1] = array1_store[indices[i]];
       array2[i+1] = array2_store[indices[i]];
    }

	delete[] indices;
	delete[] array1_store;
	delete[] array2_store;
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
