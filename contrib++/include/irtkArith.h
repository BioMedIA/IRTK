/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef IRTKARITH_H

#define IRTKARITH_H

#include <irtkImage.h>

class irtkArith
{
public:
  static int IAbs (int iValue);
  static int ICeil (float fValue);
  static int IFloor (float fValue);
  static int ISign (int iValue);

  static double Abs (double fValue);
  static double ACos (double fValue);
  static double ASin (double fValue);
  static double ATan (double fValue);
  static double ATan2 (double fY, double fX);
  static double Ceil (double fValue);
  static double Cos (double fValue);
  static double Exp (double fValue);
  static double Floor (double fValue);
  static double Log (double fValue);
  static double Pow (double kBase, double kExponent);
  static double Sign (double fValue);
  static double Sin (double fValue);
  static double Sqr (double fValue);
  static double Sqrt (double fValue);
  static double UnitRandom ();  // in [0,1]
  static double SymmetricRandom ();  // in [-1,1]

};

#endif
