/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkArith.h>

//----------------------------------------------------------------------------
int irtkArith::IAbs (int iValue)
{
  return ( iValue >= 0 ? iValue : -iValue );
}
//----------------------------------------------------------------------------
int irtkArith::ICeil (float fValue)
{
  return int(ceil(fValue));
}
//----------------------------------------------------------------------------
int irtkArith::IFloor (float fValue)
{
  return int(floor(fValue));
}
//----------------------------------------------------------------------------
int irtkArith::ISign (int iValue)
{
  return ( iValue > 0 ? +1 : ( iValue < 0 ? -1 : 0 ) );
}
//----------------------------------------------------------------------------
double irtkArith::Abs (double fValue)
{
  return double(fabs(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::ACos (double fValue)
{
  if ( -1.0 < fValue ) {
    if ( fValue < 1.0 )
      return double(acos(fValue));
    else
      return 0.0;
  } else {
    return M_PI;
  }
}
//----------------------------------------------------------------------------
double irtkArith::ASin (double fValue)
{
  if ( -1.0 < fValue ) {
    if ( fValue < 1.0 )
      return double(asin(fValue));
    else
      return -M_PI/2.0;
  } else {
    return M_PI/2.0;
  }
}
//----------------------------------------------------------------------------
double irtkArith::ATan (double fValue)
{
  return double(atan(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::ATan2 (double fY, double fX)
{
  return double(atan2(fY,fX));
}
//----------------------------------------------------------------------------
double irtkArith::Ceil (double fValue)
{
  return double(ceil(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::Cos (double fValue)
{
  return double(cos(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::Exp (double fValue)
{
  return double(exp(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::Floor (double fValue)
{
  return double(floor(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::Log (double fValue)
{
  return double(log(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::Pow (double fBase, double fExponent)
{
  return double(pow(fBase,fExponent));
}
//----------------------------------------------------------------------------
double irtkArith::Sign (double fValue)
{
  if ( fValue > 0.0 )
    return 1.0;

  if ( fValue < 0.0 )
    return -1.0;

  return 0.0;
}
//----------------------------------------------------------------------------
double irtkArith::Sin (double fValue)
{
  return double(sin(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::Sqr (double fValue)
{
  return fValue*fValue;
}
//----------------------------------------------------------------------------
double irtkArith::Sqrt (double fValue)
{
  return double(sqrt(fValue));
}
//----------------------------------------------------------------------------
double irtkArith::UnitRandom ()
{
  return double(rand())/double(RAND_MAX);
}
//----------------------------------------------------------------------------
double irtkArith::SymmetricRandom ()
{
  return 2.0*double(rand())/double(RAND_MAX) - 1.0;
}
//----------------------------------------------------------------------------
