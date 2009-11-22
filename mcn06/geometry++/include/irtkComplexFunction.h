/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCOMPLEXFUNCTION_H

#define _IRTKCOMPLEXFUNCTION_H

/**

   Complex function class.

*/

class irtkComplexFunction : public irtkObject
{

public:

  /// Virtual destructor
  virtual ~irtkComplexFunction();

  /// Pure virtual evaluation function
  virtual complex<float> Evaluate(float x, float y, float z) = 0;

};

#endif
