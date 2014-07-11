#include "gtest/gtest.h"

#include "geometry++/src/irtkMatrix.cc"

const double EPSILON = 0.0001;

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
    for (int i = 0; i < matrix.Rows(); i++) {
      for (int j = 0; j < matrix.Cols(); j++) {
	  ASSERT_EQ(expected(i, j), round(matrix(i, j)));
      }
    }
}

TEST(Geometry_irtkMatrix, SVD) {
   irtkMatrix matrix(4, 2);
   matrix.Put(0, 0, 1);
   matrix.Put(0, 1, 2);
   matrix.Put(1, 0, 3);
   matrix.Put(1, 1, 4);
   matrix.Put(2, 0, 5);
   matrix.Put(2, 1, 6);
   matrix.Put(3, 0, 7);
   matrix.Put(3, 1, 8);

   irtkMatrix u;
   irtkMatrix v;
   irtkVector w;

   matrix.SVD(u, w, v);

   irtkMatrix s(2, 2);
   s(0, 0) = w(0);
   s(0, 1) = 0;
   s(1, 0) = 0;
   s(1, 1) = w(1);

   v.Transpose();   
   
   irtkMatrix result;
   result = u * s * v;
   
   for (int i = 0; i < matrix.Rows(); i++) {
      for (int j = 0; j < matrix.Cols(); j++) {
         ASSERT_EQ(matrix(i, j), round(result(i, j)));
      }
   } 
}

TEST(Geometry_irtkMatrix, LeastSquaresFit) {
   irtkMatrix matrix(5, 2);
   matrix.Put(0, 0, 1);
   matrix.Put(0, 1, 1);
   matrix.Put(1, 0, 1);
   matrix.Put(1, 1, 2);
   matrix.Put(2, 0, 1);
   matrix.Put(2, 1, 3);
   matrix.Put(3, 0, 1);
   matrix.Put(3, 1, 4);
   matrix.Put(4, 0, 1);
   matrix.Put(4, 1, 5);

   irtkVector y(5);
   y(0) = 3;
   y(1) = 4;
   y(2) = 5;
   y(3) = 6;
   y(4) = 7; 

   irtkVector x(2);

   matrix.LeastSquaresFit(y, x);

   double expected_x[2] = {2, 1};

   for (int i = 0; i < x.Rows(); i++) {
      ASSERT_NEAR(expected_x[i], x(i), EPSILON);
   }
}

TEST(Geometry_irtkMatrix, Eigenvalues_non_symmetric_matrix) {
   irtkMatrix matrix(3, 3);
   matrix.Put(0, 0, 67);
   matrix.Put(0, 1, 120);
   matrix.Put(0, 2, 170);
   matrix.Put(1, 0, -24);
   matrix.Put(1, 1, -35);
   matrix.Put(1, 2, -40);
   matrix.Put(2, 0, -8);
   matrix.Put(2, 1, -20);
   matrix.Put(2, 2, -35);

   irtkMatrix E1, E2;
   irtkVector e;

   matrix.Eigenvalues(E1, e, E2);
    
   irtkVector exp(3);
   exp(0) = 52768.81418;
   exp(1) = 110.18485;
   exp(2) = 0.00096;

   for (int i = 0; i < e.Rows(); i++) {
      ASSERT_NEAR(exp(i), e(i), EPSILON);
   }
}

TEST(Geometry_irtkMatrix, Eigenvalues_symmetric_matrix) {
   irtkMatrix matrix(3, 3);
   matrix.Put(0, 0, 1);
   matrix.Put(0, 1, 7);
   matrix.Put(0, 2, 3);
   matrix.Put(1, 0, 7);
   matrix.Put(1, 1, 4);
   matrix.Put(1, 2, -5);
   matrix.Put(2, 0, 3);
   matrix.Put(2, 1, -5);
   matrix.Put(2, 2, 6);

   irtkMatrix E1, E2;
   irtkVector e;

   matrix.Eigenvalues(E1, e, E2);
    
   irtkVector exp(3);
   exp(0) = -7.007928;
   exp(1) = 7.03595;
   exp(2) = 10.971982;

   for (int i = 0; i < e.Rows(); i++) {
      ASSERT_NEAR(exp(i), e(i), EPSILON);
   }
}

TEST(Geometry_irtkMatrix, Eigenvalues_diag_matrix) {
   irtkMatrix matrix(3, 3);
   matrix.Put(0, 0, 1);
   matrix.Put(0, 1, 1);
   matrix.Put(0, 2, 0);
   matrix.Put(1, 0, 0);
   matrix.Put(1, 1, 1);
   matrix.Put(1, 2, 1);
   matrix.Put(2, 0, 1);
   matrix.Put(2, 1, 0);
   matrix.Put(2, 2, 1);

   irtkMatrix E1, E2;
   irtkVector e;

   ASSERT_THROW(matrix.Eigenvalues(E1, e, E2), irtkException);
}
