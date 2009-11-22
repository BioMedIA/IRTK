/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkVector3.h>

const irtkVector3 irtkVector3::ZERO(0,0,0);
const irtkVector3 irtkVector3::UNIT_X(1,0,0);
const irtkVector3 irtkVector3::UNIT_Y(0,1,0);
const irtkVector3 irtkVector3::UNIT_Z(0,0,1);

//----------------------------------------------------------------------------
irtkVector3::irtkVector3 ()
{
  // For efficiency in construction of large arrays of vectors, the
  // default constructor does not initialize the vector.
}
//----------------------------------------------------------------------------
irtkVector3::irtkVector3 (double fX, double fY, double fZ)
{
  x = fX;
  y = fY;
  z = fZ;
}
//----------------------------------------------------------------------------
irtkVector3::irtkVector3 (double afCoordinate[3])
{
  x = afCoordinate[0];
  y = afCoordinate[1];
  z = afCoordinate[2];
}
//----------------------------------------------------------------------------
irtkVector3::irtkVector3 (const irtkVector3& rkVector)
{
  x = rkVector.x;
  y = rkVector.y;
  z = rkVector.z;
}
//----------------------------------------------------------------------------
double& irtkVector3::operator[] (int i) const
{
  // assert:  0 <= i < 2; x, y, and z are packed into 3*sizeof(double)
  //          bytes
  return (double&) *(&x + i);
}
//----------------------------------------------------------------------------
irtkVector3::operator double* ()
{
  return &x;
}
//----------------------------------------------------------------------------
irtkVector3& irtkVector3::operator= (const irtkVector3& rkVector)
{
  x = rkVector.x;
  y = rkVector.y;
  z = rkVector.z;
  return *this;
}
//----------------------------------------------------------------------------
bool irtkVector3::operator== (const irtkVector3& rkVector) const
{
  return ( x == rkVector.x && y == rkVector.y && z == rkVector.z );
}
//----------------------------------------------------------------------------
bool irtkVector3::operator!= (const irtkVector3& rkVector) const
{
  return ( x != rkVector.x || y != rkVector.y || z != rkVector.z );
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::operator+ (const irtkVector3& rkVector) const
{
  irtkVector3 kSum;
  kSum.x = x + rkVector.x;
  kSum.y = y + rkVector.y;
  kSum.z = z + rkVector.z;
  return kSum;
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::operator- (const irtkVector3& rkVector) const
{
  irtkVector3 kDiff;
  kDiff.x = x - rkVector.x;
  kDiff.y = y - rkVector.y;
  kDiff.z = z - rkVector.z;
  return kDiff;
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::operator* (double fScalar) const
{
  irtkVector3 kProd;
  kProd.x = fScalar*x;
  kProd.y = fScalar*y;
  kProd.z = fScalar*z;
  return kProd;
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::operator- () const
{
  irtkVector3 kNeg;
  kNeg.x = -x;
  kNeg.y = -y;
  kNeg.z = -z;
  return kNeg;
}
//----------------------------------------------------------------------------
irtkVector3 operator* (double fScalar, const irtkVector3& rkVector)
{
  irtkVector3 kProd;
  kProd.x = fScalar*rkVector.x;
  kProd.y = fScalar*rkVector.y;
  kProd.z = fScalar*rkVector.z;
  return kProd;
}
//----------------------------------------------------------------------------
irtkVector3& irtkVector3::operator+= (const irtkVector3& rkVector)
{
  x += rkVector.x;
  y += rkVector.y;
  z += rkVector.z;
  return *this;
}
//----------------------------------------------------------------------------
irtkVector3& irtkVector3::operator-= (const irtkVector3& rkVector)
{
  x -= rkVector.x;
  y -= rkVector.y;
  z -= rkVector.z;
  return *this;
}
//----------------------------------------------------------------------------
irtkVector3& irtkVector3::operator*= (double fScalar)
{
  x *= fScalar;
  y *= fScalar;
  z *= fScalar;
  return *this;
}
//----------------------------------------------------------------------------
double irtkVector3::SquaredLength () const
{
  return x*x + y*y + z*z;
}
//----------------------------------------------------------------------------
double irtkVector3::Dot (const irtkVector3& rkVector) const
{
  return x*rkVector.x + y*rkVector.y + z*rkVector.z;
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::operator/ (double fScalar) const
{
  irtkVector3 kQuot;

  if ( fScalar != 0.0 ) {
    double fInvScalar = 1.0/fScalar;
    kQuot.x = fInvScalar*x;
    kQuot.y = fInvScalar*y;
    kQuot.z = fInvScalar*z;
    return kQuot;
  } else {
    cerr << "irtkVector3::operator/: division by 0" << endl;
    exit(1);
  }
}
//----------------------------------------------------------------------------
irtkVector3& irtkVector3::operator/= (double fScalar)
{
  if ( fScalar != 0.0 ) {
    double fInvScalar = 1.0/fScalar;
    x *= fInvScalar;
    y *= fInvScalar;
    z *= fInvScalar;
  } else {
    cerr << "irtkVector3::operator/=: division by 0" << endl;
    exit(1);
  }

  return *this;
}
//----------------------------------------------------------------------------
double irtkVector3::Length () const
{
  return irtkArith::Sqrt(x*x + y*y + z*z);
}
//----------------------------------------------------------------------------
double irtkVector3::Unitize (double fTolerance)
{
  double fLength = Length();

  if ( fLength > fTolerance ) {
    double fInvLength = 1.0/fLength;
    x *= fInvLength;
    y *= fInvLength;
    z *= fInvLength;
  } else {
    fLength = 0.0;
  }

  return fLength;
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::Cross (const irtkVector3& rkVector) const
{
  irtkVector3 kCross;

  kCross.x = y*rkVector.z - z*rkVector.y;
  kCross.y = z*rkVector.x - x*rkVector.z;
  kCross.z = x*rkVector.y - y*rkVector.x;

  return kCross;
}
//----------------------------------------------------------------------------
irtkVector3 irtkVector3::UnitCross (const irtkVector3& rkVector) const
{
  irtkVector3 kCross;

  kCross.x = y*rkVector.z - z*rkVector.y;
  kCross.y = z*rkVector.x - x*rkVector.z;
  kCross.z = x*rkVector.y - y*rkVector.x;
  kCross.Unitize();

  return kCross;
}
//----------------------------------------------------------------------------
void irtkVector3::Orthonormalize (irtkVector3 akVector[3])
{
  // If the input vectors are v0, v1, and v2, then the Gram-Schmidt
  // orthonormalization produces vectors u0, u1, and u2 as follows,
  //
  //   u0 = v0/|v0|
  //   u1 = (v1-(u0*v1)u0)/|v1-(u0*v1)u0|
  //   u2 = (v2-(u0*v2)u0-(u1*v2)u1)/|v2-(u0*v2)u0-(u1*v2)u1|
  //
  // where |A| indicates length of vector A and A*B indicates dot
  // product of vectors A and B.

  // compute u0
  akVector[0].Unitize();

  // compute u1
  double fDot0 = akVector[0].Dot(akVector[1]);
  akVector[1] -= fDot0*akVector[0];
  akVector[1].Unitize();

  // compute u2
  double fDot1 = akVector[1].Dot(akVector[2]);
  fDot0 = akVector[0].Dot(akVector[2]);
  akVector[2] -= fDot0*akVector[0] + fDot1*akVector[1];
  akVector[2].Unitize();
}
//----------------------------------------------------------------------------
void irtkVector3::GenerateOrthonormalBasis (irtkVector3& rkU, irtkVector3& rkV,
    irtkVector3& rkW, bool bUnitLengthW)
{
  if ( !bUnitLengthW )
    rkW.Unitize();

  if ( irtkArith::Abs(rkW.x) >= irtkArith::Abs(rkW.y)
       &&   irtkArith::Abs(rkW.x) >= irtkArith::Abs(rkW.z) ) {
    rkU.x = -rkW.y;
    rkU.y = +rkW.x;
    rkU.z = 0.0;
  } else {
    rkU.x = 0.0;
    rkU.y = +rkW.z;
    rkU.z = -rkW.y;
  }

  rkU.Unitize();
  rkV = rkW.Cross(rkU);
}
//----------------------------------------------------------------------------
