/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageRegistration.h 307 2011-04-10 22:39:07Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2011-04-10 23:39:07 +0100 (Sun, 10 Apr 2011) $
  Version   : $Revision: 307 $
  Changes   : $Author: ws207 $

=========================================================================*/

#ifndef _IRTKIMAGEREGISTRATIONWITHPADDING_H

#define _IRTKIMAGEREGISTRATIONWITHPADDING_H

/**
 * Generic for image registration extended by source padding
**/

class irtkImageRegistrationWithPadding : public irtkImageRegistration
{

protected:

  /// Padding value of source image
  short  _SourcePadding;
  
  //irtkGreyImage *tmp_target, *tmp_source;

  /// Overload initial set up for the registration at a multiresolution level
  virtual void Initialize(int);

public:
  irtkImageRegistrationWithPadding();
};

#include <irtkImageRigidRegistrationWithPadding.h>

#endif
