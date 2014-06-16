#include "gtest/gtest.h"

#include "packages/segmentation/src/irtkGraphCutSegmentation_4D.cc"

TEST(Packages_Segmentation_irtkGraphCutSegmentation_4D, sort_null_array) {
  float *array = NULL;

  sort(array, array);

  ASSERT_EQ(0, array);
}

TEST(Packages_Segmentation_irtkGraphCutSegmentation_4D, sort_single_entry_array) {
  float array[1] = {5};

  sort(array, array+1);

  ASSERT_EQ(5, array[0]);
}

TEST(Packages_Segmentation_irtkGraphCutSegmentation_4D, sort_simple_array) {
  float array[5] = {1.123, 6.2134, 2.234, 512.2, 3.123};
  float sorted_array[5] = {1.123, 2.234, 3.123, 6.2134, 512.2};

  sort(array, array+5);

  for (int i = 0; i < 5; i++)
    ASSERT_EQ(sorted_array[i], array[i]); 	
}

TEST(Packages_Segmentation_irtkGraphCutSegmentation_4D, sort_array_same_elements) {
  float array[5] = {1, 10, 1, 2, 2};
  float sorted_array[5] = {1, 1, 2, 2, 10};

  sort(array, array+5);

  for (int i = 0; i < 5; i++)
    ASSERT_EQ(sorted_array[i], array[i]); 
}

