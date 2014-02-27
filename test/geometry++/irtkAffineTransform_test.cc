#include "gtest/gtest.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkPointSet.h>


static boost::random::mt19937 rng;
static boost::uniform_real<> uniform_50 (-50, 51);
static boost::uniform_real<> uniform_180(-180, 181);
static boost::uniform_real<> uniform_45 (-45, 46);

static const double EPSILON = 0.0000001;
static const double DELTA   = 0.001;

TEST(Geometry_irtkAffineTransform, Compute_error) {
  
  irtkPointSet pset;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > engine_50(rng, uniform_50);
  
  // Generate 100 different points
  for (int i = 0; i < 100; i++) {
    pset.Add(irtkPoint(engine_50(), engine_50(), engine_50()));
  }
  
  irtkAffineTransformation trans;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > engine_180(rng, uniform_180);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > engine_45(rng, uniform_45);
  
  // Loop over 100 random transformations
  for (int i = 0; i < 100; i++) {
    
    // Generate random transformation parameters    
    trans.Put(0,  engine_50 () );
    trans.Put(1,  engine_50 () );
    trans.Put(2,  engine_50 () );
    trans.Put(3,  engine_180() );
    trans.Put(4,  engine_180() );
    trans.Put(5,  engine_180() );
    trans.Put(6,  engine_50 () );
    trans.Put(7,  engine_50 () );
    trans.Put(8,  engine_50 () );
    trans.Put(9,  engine_45 () );
    trans.Put(10, engine_45 () );
    trans.Put(11, engine_45 () );

    // Loop over dofs
    for (int j = 0; j < trans.NumberOfDOFs(); j++) {
      
      // Loop over random points
      for (int k = 0; k < 100; k++) {
        
        double tmp;
        double jac1[3];
        double jac2[3];
        
        irtkPoint p;
        irtkPoint p1;
        irtkPoint p2;
        
        p = pset(k);
        
        // Compute Jacobian analytically
        trans.JacobianDOFs(jac1, j, p._x, p._y, p._z);
        
        // Compute Jacobian numerically
        tmp = trans.Get(j);
        trans.Put(j, tmp + DELTA);
        p1 = p;
        trans.Transform(p1._x, p1._y, p1._z);
        trans.Put(j, tmp - DELTA);
        p2 = p;
        trans.Transform(p2._x, p2._y, p2._z);
        trans.Put(j, tmp);
        jac2[0] = (p1._x - p2._x) / (2.0 * DELTA);
        jac2[1] = (p1._y - p2._y) / (2.0 * DELTA);
        jac2[2] = (p1._z - p2._z) / (2.0 * DELTA);
        
        // Compute error
        double a, b, c;
        a = sqrt((jac1[0] - jac2[0])*(jac1[0] - jac2[0]));
        b = sqrt((jac1[1] - jac2[1])*(jac1[1] - jac2[1]));
        c = sqrt((jac1[2] - jac2[2])*(jac1[2] - jac2[2]));
         
        ASSERT_FALSE (a > EPSILON);
        ASSERT_FALSE (b > EPSILON);
        ASSERT_FALSE (c > EPSILON);
      }
    }
  }
}
