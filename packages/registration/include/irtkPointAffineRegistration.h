/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPOINTAFFINEREGISTRATION_H

#define _IRTKPOINTAFFINEREGISTRATION_H

#include <irtkImage.h>

#include <irtkTransformation.h>

/**
 * Filter for point-based affine registration.
 *
 * This class implements a registration filter for point-based registration
 * of two sets of points. As an output, it returns the affine transformation
 * coefficients in a irtkAffineTransformation class.
 *
*/

class irtkPointAffineRegistration : public irtkPointRigidRegistration
{

public:

  /// Constructor
  irtkPointAffineRegistration();

  /// Destructor
  virtual ~irtkPointAffineRegistration();

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Run the filter
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass();

};

#endif



