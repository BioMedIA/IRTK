/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHOMOGENEOUSIMAGETRANSFORMATION_H

#define _IRTKHOMOGENEOUSIMAGETRANSFORMATION_H

#include <irtkImage.h>

#include <irtkTransformation.h>

class irtkImageHomogeneousTransformation  : public irtkImageTransformation
{

public:

  /** Constructor. This constructs an transformation filter with a given
   *  interpolation mode and padding value. By default the interpolation
   *  mode is set to trilinear.
   */
  irtkImageHomogeneousTransformation();

  /// Destructor
  virtual ~irtkImageHomogeneousTransformation();

  /// Sets transformation
  virtual void SetTransformation(irtkTransformation *);

  /// Runs the filter
  virtual void Run();

};

#endif
