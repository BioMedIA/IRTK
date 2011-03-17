/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkRician.h 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#ifndef _IRTKRician_H

#define _IRTKRician_H

#include <irtkImage.h>

#include <nr.h>

/**

multivariate gaussian probability distribution

*/

class irtkRician : public irtkObject
{

protected:

  double _mi;
  double _sigma;
  double _norm;

  double Bessel(const double &);

public:

  void Initialise(const double &mi, const double &sigma);

  /// Approximate parameters from mean and std given by initialise
  /// if initialise with rician parameters not mean and std then this function should not be called
  void Approximate();

  double Evaluate(const double &x);

  GetMacro(norm,double);
};

#endif
