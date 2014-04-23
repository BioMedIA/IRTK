#include "../include/makevolume_core.h"

#include <irtkImage.h>

#include <algorithm>

long unsigned* indexing(long unsigned n, float* array) {
  // initialize vector with the array values
  // pointers can be used like any other iterators
  vector<float> vec(array, array+n);

  // initialize indices vector
  vector<long unsigned> idx;
  idx.resize(n);
  for (long unsigned i = 0; i < n; i++) idx[i] = i;

  // do the sorting
  sort(idx.begin(), idx.end(), index_sorting<vector<float> &>(vec));

  // fill the array to maintain function signature
  long unsigned* index = new long unsigned[n];
  for (long unsigned i = 0; i < n; i++) index[i] = idx[i]+1; 

  return index;
}
