/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#ifndef _IRTKMULTIPLESURFACEAFFINEREGISTRATION_H

#define _IRTKMULTIPLESURFACEAFFINEREGISTRATION_H

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkPointRegistration.h>
#include <irtkMultipleSurfaceRegistration.h>

#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <irtkLocator.h>
#include <vtkFeatureEdges.h>

/**
 * Filter for affine registration based on surfaces.
 *
 */

class irtkMultipleSurfaceAffineRegistration : public irtkMultipleSurfaceRegistration
{

protected:

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

public:

  /// Constructor
  irtkMultipleSurfaceAffineRegistration();

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

#endif

#endif
