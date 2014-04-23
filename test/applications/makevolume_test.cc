#include "gtest/gtest.h"
#include "applications/include/makevolume_core.h"

TEST(Applications_Makevolume, indexing_single_entry_array) {
  long unsigned index[1] = {1};
  float array[100] = {5};

  ASSERT_EQ(index[0], indexing(1, array)[0]);
}

TEST(Applications_Makevolume, indexing_simple_array) {
  long unsigned index[5] = {1, 2, 3, 4, 5};
  float array[5] = {5, 6, 7, 8, 9};

  long unsigned* result = indexing(5, array);
  
  for (int i = 0; i < 5; i++)
    ASSERT_EQ(index[i], result[i]);
}

TEST(Applications_Makevolume, indexing_reverse_array) {
  long unsigned index[4] = {4, 3, 2, 1};
  float array[4] = {1000.97, 234.1, 200, 1.32435};

  long unsigned* result = indexing(4, array);

  for (int i = 0; i < 4; i++)
    ASSERT_EQ(index[i], result[i]);
}

TEST(Applications_Makevolume, indexing_random_array) {
  long unsigned index[10] = {1, 3, 6, 4, 7, 2, 8, 5, 10, 9};
  float array[10] = {3.34, 8.234, 4.1, 6.1324, 10.354235, 5.12, 7.123, 9.2345, 12.32, 11.0};

  long unsigned* result = indexing(10, array);

  for (int i = 0; i < 10; i++)
    ASSERT_EQ(index[i], result[i]);
} 
