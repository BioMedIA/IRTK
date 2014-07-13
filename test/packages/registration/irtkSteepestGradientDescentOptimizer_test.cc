#include "gtest/gtest.h"

#include <irtkImage.h>
#include <irtkRegistration.h>
#include <irtkTransformation.h>

const double EPSILON = 0.0001;

TEST(Packages_Registration_irtkSteepestGradientDescientOptimizer, Run) {
   irtkGreyImage *target = new irtkGreyImage;
   irtkGreyImage *source = new irtkGreyImage;
   target->Read("/vol/vipdata/data/brain/adni/images/ADNIGO/ADNI_002_S_0685_MR_MT1__GradWarp__N3m_Br_20120322163357822_S89145_I291874.nii.gz");
   source->Read("/vol/vipdata/data/brain/adni/images/ADNIGO/ADNI_003_S_0907_MR_MT1__GradWarp__N3m_Br_20120322164718634_S99132_I291887.nii.gz");

   irtkRigidTransformation transformation;

   irtkImageRigidRegistration registration;
   registration.SetInput(target, source);
   registration.SetOutput(&transformation);
   registration.SetOptimizationMethod(SteepestGradientDescent);
   registration.Run();

   for (int i = 0; i < transformation.NumberOfDOFs(); i++) {
      ASSERT_NEAR(0.0, transformation.Get(i), EPSILON);
   }

   delete target;
   delete source;
}
