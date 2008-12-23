/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPOINTRIGIDREGISTRATION_H

#define _IRTKPOINTRIGIDREGISTRATION_H

#include <irtkImage.h>

#include <irtkTransformation.h>

/**
 * Filter for point-based registration.
 *
 * This class implements a registration filter for point-based registration
 * of two sets of points. As an output, it returns the rigid transformation
 * coefficients in a irtkRigidTransformation class.
 *
*/

class irtkPointRigidRegistration : public irtkPointRegistration
{

protected:

  /// Initial set up for the registration
  virtual void Initialize();

  /// Optimize registration using closed form solution
  virtual void ClosedFormOptimizer();

public:

  /// Constructor
  irtkPointRigidRegistration();

  /// Destructor
  virtual ~irtkPointRigidRegistration();

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Run the filter
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass();

};

#endif
