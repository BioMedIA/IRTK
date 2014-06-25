#include "gtest/gtest.h"
#include "packages/segmentation/include/irtkRician.h"

static const double EPSILON = 0.00001;

TEST(Packages_Segmentation_irtkRician, Approximate_0_1) {
   irtkRician *rician = new irtkRician;
   
   rician->Initialise(0, 1);
   rician->Approximate();
  
   double norm;
   norm = rician->Getnorm();
   
   ASSERT_EQ(0.429, norm);

   delete rician;
}

TEST(Packages_Segmentation_irtkRician, Approximate_5_5) {
   irtkRician *rician = new irtkRician;

   rician->Initialise(5, 5);
   rician->Approximate();

   double norm, expected_norm = 0.151337;
   norm = rician->Getnorm();

   ASSERT_NEAR(expected_norm, norm, EPSILON);

   delete rician;
}

TEST(Packages_Segmentation_irtkRician, Approximate_5_10) {
   irtkRician *rician = new irtkRician;

   rician->Initialise(5, 10);
   rician->Approximate();

   double norm;
   norm = rician->Getnorm();

   ASSERT_EQ(0.0429, norm);

   delete rician;
}
