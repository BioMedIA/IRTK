#include "../include/makevolume_core.h"

#include <irtkImage.h>
#include <nr.h>

long unsigned* indexing(unsigned long n, float* array) {
  long unsigned* index = new unsigned long[n];

  indexx(n, array-1, index-1);

  return index;
}
