#include "gtest/gtest.h"

#include <sys/types.h>
#include <sys/time.h>

#include <cmath>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


static boost::mt19937 rng;
static boost::normal_distribution<> nd(0, 1);

static const double EPSILON = 0.001;
static const int SIZE = 1000000;

TEST(Image_irtkGaussianNoise, Mean) {
   double expected_mean = 0.0, mean = 0.0;

   timeval tv;
   gettimeofday(&tv, NULL);
   long init = tv.tv_usec;

   rng.seed(init);
   boost::variate_generator<boost::mt19937&, 
                            boost::normal_distribution<> > var_nor(rng, nd);

   for (int i = 0; i < SIZE; i++) {
      mean += var_nor();
   }
   mean /= SIZE;

   ASSERT_NEAR(expected_mean, mean, EPSILON);
}

TEST(Image_irtkGaussianNoise, StandardDeviation) {
   double mean = 0.0, expected_std = 1.0, std = 0.0;
   double nums[SIZE];

   timeval tv;
   gettimeofday(&tv, NULL);
   long init = tv.tv_sec; 

   rng.seed(init);
   boost::variate_generator<boost::mt19937&,
			    boost::normal_distribution<> > var_nor(rng, nd);

   for (int i = 0; i < SIZE; i++) {
       nums[i] = var_nor();
       mean += nums[i];
   }

   mean /= SIZE;

   for (int i = 0; i < SIZE; i++) {
       std += pow(nums[i]-mean, 2);
   }

   std = sqrt(std / SIZE);

   ASSERT_NEAR(expected_std, std, EPSILON);
}
