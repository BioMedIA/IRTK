#include "gtest/gtest.h"

#include "geometry++/include/irtkPointSet.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


static boost::random::mt19937 rng;
static boost::uniform_real<> uniform_50(-50, 51);

static const double EPSILON = 0.0001;

TEST(Geometry_irtkPointSet, Center_of_gravity) {

  
  const irtkPoint expectedCentreOfGravity(0.0961255, 1.52611, 5.24152);
  irtkPointSet pset;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > engine(rng, uniform_50);
  
  // Generate 100 different points
  for (int i = 0; i < 100; i++) {
    pset.Add(irtkPoint(engine(), engine(), engine()));
  }
  
  ASSERT_NEAR(expectedCentreOfGravity._x, pset.CenterOfGravity()._x, EPSILON);
  ASSERT_NEAR(expectedCentreOfGravity._y, pset.CenterOfGravity()._y, EPSILON);
  ASSERT_NEAR(expectedCentreOfGravity._z, pset.CenterOfGravity()._z, EPSILON);
}

TEST(Geometry_irtkPointSet, Standard_deviation_ellipsoid) {
  const irtkPoint expectedStandardDeviationEllipsoid (28.0534, 29.1969, 26.7527);
  
  irtkPointSet pset;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > engine(rng, uniform_50);
  
  // Generate 100 different points
  for (int i = 0; i < 100; i++) {
    pset.Add(irtkPoint(engine(), engine(), engine()));
  }
  
  ASSERT_NEAR(expectedStandardDeviationEllipsoid._x, pset.StandardDeviationEllipsoid()._x, EPSILON);
  ASSERT_NEAR(expectedStandardDeviationEllipsoid._y, pset.StandardDeviationEllipsoid()._y, EPSILON);
  ASSERT_NEAR(expectedStandardDeviationEllipsoid._z, pset.StandardDeviationEllipsoid()._z, EPSILON);
}
