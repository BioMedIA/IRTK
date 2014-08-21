#include "gtest/gtest.h"

#include <irtkImage.h>
#include <irtkRegistration.h>
#include <irtkTransformation.h>

#ifdef _WIN32
#include <Windows.h>
#endif

const double EPSILON = 0.0001;

TEST(Packages_Registration_irtkConjugateGradientDescentOptimizer, Run) {
   irtkGreyImage *target = new irtkGreyImage;
   irtkGreyImage *source = new irtkGreyImage;

   std::string image_path = "";

#ifdef __linux__
   const int buffer_size = 1024;
   char buffer[buffer_size];
   ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer)-1);

   std::string key ("irtk/");

   if (len != -1) {
       buffer[len] = '\0';
       std::string app_path = std::string(buffer);
       unsigned found = app_path.rfind(key);
       image_path = app_path.substr(0, found) + key + "test/image_files/";
   }
   else {
       std::cerr << "Could not get application path" << std::endl;
       exit(1);
   }
#elif _WIN32
   std::string key ("irtk\\");

    char buffer[MAX_PATH];
    GetModuleFileName( NULL, buffer, MAX_PATH );
    std::string app_path = std::string(buffer);
    unsigned found = app_path.rfind(key);
    root_path = app_path.substr(0, found) + key + "test\\image_files\\";
#endif

   target->Read((image_path + "target.nii.gz").c_str());
   source->Read((image_path + "source.nii.gz").c_str());

   irtkRigidTransformation transformation;

   irtkImageRigidRegistration registration;
   registration.SetInput(target, source);
   registration.SetOutput(&transformation);
   registration.SetOptimizationMethod(ConjugateGradientDescent);
   registration.Run(); 

   for (int i = 0; i < transformation.NumberOfDOFs(); i++) {
      ASSERT_NEAR(0, transformation.Get(i), EPSILON);
   }
   //transformation.Print();
   
   delete target;
   delete source;
}
