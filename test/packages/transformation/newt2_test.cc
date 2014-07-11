#include "gtest/gtest.h"

//#include <nr.h>
//#include <nrutil.h>
#include "packages/transformation/src/newt2.cc"
#include <irtkTransformation.h>

const double EPSILON = 0.001;

TEST(Packages_Transformation_newt2, irtkMultiLevelFreeFormTransformation_Inverse)
{
   irtkCifstream iStream;
   iStream.Open("/vol/medic02/users/sp2010/PhD/Tumor/dofs/gradientnreg2-NGP-HG01.dof.gz");
   irtkMultiLevelFreeFormTransformation transformation;
   transformation.Read(iStream);

   double x = 10, y = 10, z = 10;
   transformation.Inverse(x, y, z);

   ASSERT_NEAR(1.16573, x, EPSILON);
   ASSERT_NEAR(23.9054, y, EPSILON);
   ASSERT_NEAR(-17.255, z, EPSILON);
}

