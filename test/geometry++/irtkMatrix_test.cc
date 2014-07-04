#include "gtest/gtest.h"

#include "geometry++/src/irtkMatrix.cc"

TEST(Geometry_irtkMatrix, Det_3x3) {
    irtkMatrix matrix(3, 3);

    // First row
    matrix.Put(0, 0, 6);
    matrix.Put(0, 1, 1);
    matrix.Put(0, 2, 1);

    // Second row
    matrix.Put(1, 0, 4);
    matrix.Put(1, 1, -2);
    matrix.Put(1, 2, 5);

    // Third row
    matrix.Put(2, 0, 2);
    matrix.Put(2, 1, 8);
    matrix.Put(2, 2, 7);

    ASSERT_EQ(-306, matrix.Det());
}

TEST(Geometry_irtkMatrix, Invert_3x3) {
    irtkMatrix matrix(3, 3);
    matrix.Put(0, 0, 1);
    matrix.Put(0, 1, 2);
    matrix.Put(0, 2, 3);
    matrix.Put(1, 0, 0);
    matrix.Put(1, 1, 1);
    matrix.Put(1, 2, 4);
    matrix.Put(2, 0, 5);
    matrix.Put(2, 1, 6);
    matrix.Put(2, 2, 0);

    irtkMatrix expected(3, 3);
    expected.Put(0, 0,-24);
    expected.Put(0, 1, 18);
    expected.Put(0, 2,  5);
    expected.Put(1, 0, 20);
    expected.Put(1, 1,-15);
    expected.Put(1, 2, -4);
    expected.Put(2, 0, -5);
    expected.Put(2, 1,  4);
    expected.Put(2, 2,  1);

    matrix.Invert();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
	  ASSERT_EQ(expected(i, j), round(matrix(i, j)));
      }
    }
}

TEST(Geometry_irtkMatrix, Adjugate_3x3) {
    irtkMatrix matrix(3, 3);
    matrix.Put(0, 0, -3);
    matrix.Put(0, 1,  2);
    matrix.Put(0, 2, -5);
    matrix.Put(1, 0, -1);
    matrix.Put(1, 1,  0);
    matrix.Put(1, 2, -2);
    matrix.Put(2, 0,  3);
    matrix.Put(2, 1, -4);
    matrix.Put(2, 2,  1);

    irtkMatrix expected(3, 3);
    expected.Put(0, 0, -8);
    expected.Put(0, 1, 18);
    expected.Put(0, 2, -4);
    expected.Put(1, 0, -5);
    expected.Put(1, 1, 12);
    expected.Put(1, 2, -1);
    expected.Put(2, 0,  4);
    expected.Put(2, 1, -6);
    expected.Put(2, 2,  2);

    double d; 
    matrix.Adjugate(d);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
	  ASSERT_EQ(expected(i, j), round(matrix(i, j)));
      }
    }
}

TEST(Geometry_irtkMatrix, SVD) {
   irtkMatrix matrix(4, 2);
   matrix.Put(0, 0, 2);
   matrix.Put(0, 1, 4);
   matrix.Put(1, 0, 1);
   matrix.Put(1, 1, 3);
   matrix.Put(2, 0, 0);
   matrix.Put(2, 1, 0);
   matrix.Put(3, 0, 0);
   matrix.Put(3, 1, 0);  
}
