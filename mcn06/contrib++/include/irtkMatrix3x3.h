/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef IRTKMATRIX3X3_H

#define IRTKMATRIX3X3_H

#include <irtkVector3.h>

// NOTE.  The (x,y,z) coordinate system is assumed to be right-handed.
// Coordinate axis rotation matrices are of the form
//   RX =    1       0       0
//           0     cos(t) -sin(t)
//           0     sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the yz-plane
//   RY =  cos(t)    0     sin(t)
//           0       1       0
//        -sin(t)    0     cos(t)
// where t > 0 indicates a counterclockwise rotation in the zx-plane
//   RZ =  cos(t) -sin(t)    0
//         sin(t)  cos(t)    0
//           0       0       1
// where t > 0 indicates a counterclockwise rotation in the xy-plane.


class irtkMatrix3x3
{
public:
  // construction
  irtkMatrix3x3 ();
  irtkMatrix3x3 (const double aafEntry[3][3]);
  irtkMatrix3x3 (const irtkMatrix3x3& rkMatrix);
  irtkMatrix3x3 (double fEntry00, double fEntry01, double fEntry02,
                 double fEntry10, double fEntry11, double fEntry12,
                 double fEntry20, double fEntry21, double fEntry22);

  // member access, allows use of construct mat[r][c]
  double* operator[] (int iRow) const;
  operator double* ();
  irtkVector3 GetColumn (int iCol) const;

  // assignment and comparison
  irtkMatrix3x3& operator= (const irtkMatrix3x3& rkMatrix);
  bool operator== (const irtkMatrix3x3& rkMatrix) const;
  bool operator!= (const irtkMatrix3x3& rkMatrix) const;

  // arithmetic operations
  irtkMatrix3x3 operator+ (const irtkMatrix3x3& rkMatrix) const;
  irtkMatrix3x3 operator- (const irtkMatrix3x3& rkMatrix) const;
  irtkMatrix3x3 operator* (const irtkMatrix3x3& rkMatrix) const;
  irtkMatrix3x3 operator- () const;

  // matrix * vector [3x3 * 3x1 = 3x1]
  irtkVector3 operator* (const irtkVector3& rkVector) const;

  // vector * matrix [1x3 * 3x3 = 1x3]
  friend irtkVector3 operator* (const irtkVector3& rkVector,
                                const irtkMatrix3x3& rkMatrix);

  // matrix * scalar
  irtkMatrix3x3 operator* (double fScalar) const;

  // scalar * matrix
  friend irtkMatrix3x3 operator* (double fScalar, const irtkMatrix3x3& rkMatrix);

  // utilities
  irtkMatrix3x3 Transpose () const;
  bool Inverse (irtkMatrix3x3& rkInverse, double fTolerance = 1e-06) const;
  irtkMatrix3x3 Inverse (double fTolerance = 1e-06) const;
  double Determinant () const;

  // singular value decomposition
  void SingularValueDecomposition (irtkMatrix3x3& rkL, irtkVector3& rkS,
                                   irtkMatrix3x3& rkR) const;
  void SingularValueComposition (const irtkMatrix3x3& rkL,
                                 const irtkVector3& rkS, const irtkMatrix3x3& rkR);

  // Gram-Schmidt orthonormalization (applied to columns of rotation matrix)
  void Orthonormalize ();

  // orthogonal Q, diagonal D, upper triangular U stored as (u01,u02,u12)
  void QDUDecomposition (irtkMatrix3x3& rkQ, irtkVector3& rkD,
                         irtkVector3& rkU) const;

  double SpectralNorm () const;

  // matrix must be orthonormal
  void ToAxisAngle (irtkVector3& rkAxis, double& rfRadians) const;
  void FromAxisAngle (const irtkVector3& rkAxis, double fRadians);

  // The matrix must be orthonormal.  The decomposition is yaw*pitch*roll
  // where yaw is rotation about the Up vector, pitch is rotation about the
  // Right axis, and roll is rotation about the Direction axis.
  bool ToEulerAnglesXYZ (float& rfYAngle, float& rfPAngle,
                         float& rfRAngle) const;
  bool ToEulerAnglesXZY (float& rfYAngle, float& rfPAngle,
                         float& rfRAngle) const;
  bool ToEulerAnglesYXZ (float& rfYAngle, float& rfPAngle,
                         float& rfRAngle) const;
  bool ToEulerAnglesYZX (float& rfYAngle, float& rfPAngle,
                         float& rfRAngle) const;
  bool ToEulerAnglesZXY (float& rfYAngle, float& rfPAngle,
                         float& rfRAngle) const;
  bool ToEulerAnglesZYX (float& rfYAngle, float& rfPAngle,
                         float& rfRAngle) const;
  void FromEulerAnglesXYZ (float fYAngle, float fPAngle, float fRAngle);
  void FromEulerAnglesXZY (float fYAngle, float fPAngle, float fRAngle);
  void FromEulerAnglesYXZ (float fYAngle, float fPAngle, float fRAngle);
  void FromEulerAnglesYZX (float fYAngle, float fPAngle, float fRAngle);
  void FromEulerAnglesZXY (float fYAngle, float fPAngle, float fRAngle);
  void FromEulerAnglesZYX (float fYAngle, float fPAngle, float fRAngle);

  // eigensolver, matrix must be symmetric
  void EigenSolveSymmetric (double afEigenvalue[3],
                            irtkVector3 akEigenvector[3]) const;

  static void TensorProduct (const irtkVector3& rkU, const irtkVector3& rkV,
                             irtkMatrix3x3& rkProduct);

  static const double EPSILON;
  static const irtkMatrix3x3 ZERO;
  static const irtkMatrix3x3 IDENTITY;

protected:
  // support for eigensolver
  void Tridiagonal (double afDiag[3], double afSubDiag[3]);
  bool QLAlgorithm (double afDiag[3], double afSubDiag[3]);

  // support for singular value decomposition
  static const double ms_fSvdEpsilon;
  static const int ms_iSvdMaxIterations;
  static void Bidiagonalize (irtkMatrix3x3& kA, irtkMatrix3x3& kL,
                             irtkMatrix3x3& kR);
  static void GolubKahanStep (irtkMatrix3x3& kA, irtkMatrix3x3& kL,
                              irtkMatrix3x3& kR);

  // support for spectral norm
  static double MaxCubicRoot (double afCoeff[3]);

  double m_aafEntry[3][3];
};

#endif
