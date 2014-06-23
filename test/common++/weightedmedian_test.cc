#include "gtest/gtest.h"

#include "common++/include/irtkCommon.h"

TEST(Common_weightedmedian, weightedmedian_single_entry_arrays) {
  int index = 2;
  double lambda = 0.5;
  double f = 2;

  float neighbors[] = {10};
  float weight[] = {0.5};

  double result = weightedmedian(index, lambda, f, neighbors, weight);

  ASSERT_EQ(6, result);
}

TEST(Common_weightedmedian, weightedmedian_simple_arrays) {
  int index = 6;
  double lambda = 1;
  double f = 1;

  float neighbors[] = {5, 3, 2, 6, 1};
  float weight[] = {-10, 8, -1, 3, 5};

  double result = weightedmedian(index, lambda, f, neighbors, weight);

  ASSERT_EQ(1, result);
}
